#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
using namespace std;

// Structure for interval information in final output
struct IntervalInfo {
    int end;
    int size;
    int count;
    bool is_rc;
};

// Structure for an interval during processing
struct Interval {
    int start;
    int end;
    bool is_rc;
};

class DNASequenceDivider {
private:
    string reference;
    string query;
    int n; // query length
    int m; // reference length
    unordered_map<string, pair<int, bool>> substring_map; // substring -> {start, is_rc}

    // Compute reverse complement of a string
    string get_reverse_complement(const string& s) {
        string res = s;
        reverse(res.begin(), res.end());
        for (char& c : res) {
            if (c == 'A') c = 'T';
            else if (c == 'T') c = 'A';
            else if (c == 'C') c = 'G';
            else if (c == 'G') c = 'C';
        }
        return res;
    }

    // Preprocess reference and its reverse complement
    void preprocess() {
        // Process reference substrings
        for (int i = 0; i < m; i++) {
            string sub = "";
            for (int j = i; j < m; j++) {
                sub += reference[j];
                substring_map[sub] = {i, false};
            }
        }

        // Process reverse complement substrings
        string reference_rc = get_reverse_complement(reference);
        for (int i = 0; i < m; i++) {
            string sub = "";
            for (int j = i; j < m; j++) {
                sub += reference_rc[j];
                int len = j - i + 1;
                int original_start = m - len - i;
                if (substring_map.find(sub) == substring_map.end()) {
                    substring_map[sub] = {original_start, true};
                }
            }
        }
    }

    // Compute dynamic programming solution
    vector<Interval> compute_dp() {
        vector<int> dp(n + 1, n + 1);
        dp[n] = 0;
        vector<pair<int, pair<int, bool>>> trace(n + 1, {-1, {-1, false}});
        unordered_map<string, int> interval_count;

        for (int i = n - 1; i >= 0; i--) {
            string sub = "";
            int min_segments = n + 1;
            pair<int, pair<int, bool>> best_trace = {-1, {-1, false}};
            for (int j = i; j < n; j++) {
                sub += query[j];
                if (substring_map.find(sub) != substring_map.end()) {
                    int next_pos = j + 1;
                    int segments = 1 + dp[next_pos];
                    int start = substring_map[sub].first;
                    bool is_rc = substring_map[sub].second;
                    int len = next_pos - i;
                    string interval_key = to_string(start) + "_" + to_string(start + len - 1);

                    if (segments < min_segments) {
                        min_segments = segments;
                        best_trace = {next_pos, {start, is_rc}};
                    } else if (segments == min_segments) {
                        int current_count = interval_count[interval_key];
                        string best_interval = to_string(best_trace.second.first) + "_" +
                                             to_string(best_trace.second.first + (next_pos - i) - 1);
                        int best_count = interval_count[best_interval];
                        if (current_count > best_count) {
                            best_trace = {next_pos, {start, is_rc}};
                        }
                    }
                }
            }
            if (min_segments <= n) {
                dp[i] = min_segments;
                trace[i] = best_trace;
                int start = best_trace.second.first;
                int len = best_trace.first - i;
                string interval_key = to_string(start) + "_" + to_string(start + len - 1);
                interval_count[interval_key]++;
            }
        }

        // Reconstruct intervals
        return reconstruct_intervals(trace);
    }

    // Reconstruct intervals from trace
    vector<Interval> reconstruct_intervals(const vector<pair<int, pair<int, bool>>>& trace) {
        vector<Interval> intervals;
        int pos = 0;
        while (pos < n) {
            int next = trace[pos].first;
            int start = trace[pos].second.first;
            bool is_rc = trace[pos].second.second;
            int len = next - pos;
            int end = start + len - 1;
            intervals.push_back({start + 1, end + 1, is_rc}); // 1-based indexing
            pos = next;
        }
        return intervals;
    }

    // Compute interval statistics
    vector<IntervalInfo> compute_interval_stats(const vector<Interval>& intervals) {
        unordered_map<string, int> occurrence_count; // key: "start_end_is_rc"
        for (const auto& interval : intervals) {
            string key = to_string(interval.start) + "_" + 
                        to_string(interval.end) + "_" + 
                        (interval.is_rc ? "1" : "0");
            occurrence_count[key]++;
        }

        vector<IntervalInfo> stats;
        unordered_map<string, IntervalInfo> interval_stats; // key: "end_size_is_rc"

        for (const auto& small : intervals) {
            int s_start = small.start;
            int s_end = small.end;
            bool s_rc = small.is_rc;
            bool is_contained = false;

            for (const auto& big : intervals) {
                int b_start = big.start;
                int b_end = big.end;
                bool b_rc = big.is_rc;
                if (b_start <= s_start && s_end <= b_end && 
                    !(b_start == s_start && b_end == s_end && b_rc == s_rc)) {
                    is_contained = true;
                    break;
                }
            }

            if (is_contained) {
                int size = s_end - s_start + 1;
                string occur_key = to_string(s_start) + "_" + to_string(s_end) + "_" + (s_rc ? "1" : "0");
                int repeat_count = occurrence_count[occur_key];
                string stat_key = to_string(s_end) + "_" + to_string(size) + "_" + (s_rc ? "1" : "0");
                if (interval_stats.find(stat_key) == interval_stats.end()) {
                    interval_stats[stat_key] = {s_end, size, repeat_count, s_rc};
                }
            }
        }

        for (const auto& [key, info] : interval_stats) {
            stats.push_back(info);
        }
        return stats;
    }

    // Print results in required format
    void print_results(const vector<IntervalInfo>& stats) {
        for (const auto& info : stats) {
            cout << "POS in REF:" << info.end 
                 << " repeat size:" << info.size 
                 << " repeat count:" << info.count 
                 << " inverse:" << (info.is_rc ? "yes" : "no") << endl;
        }
    }

public:
    DNASequenceDivider(const string& ref, const string& qry) 
        : reference(ref), query(qry), n(qry.length()), m(ref.length()) {}

    void solve() {
        preprocess();
        vector<Interval> intervals = compute_dp();
        vector<IntervalInfo> stats = compute_interval_stats(intervals);
        print_results(stats);
    }
};

int main() {
    string reference, query;
    cin >> reference >> query;
    DNASequenceDivider divider(reference, query);
    divider.solve();
    return 0;
}
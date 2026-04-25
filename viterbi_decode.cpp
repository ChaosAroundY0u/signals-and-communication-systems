#include "mex.h"
#include <algorithm>
#include <climits>
#include <vector>

class ViterbiDecoder80211a {
private:
    static constexpr int K = 7;
    static constexpr int NUM_STATES = 64;
    static constexpr int TAIL_BITS = 6;
    
    uint8_t expected_output[NUM_STATES][2];
    
    void computeExpectedOutput() {
        for (int state = 0; state < NUM_STATES; state++) {
            for (int input = 0; input < 2; input++) {
                uint8_t reg[6];
                for (int i = 0; i < 6; i++) {
                    reg[5 - i] = (state >> i) & 1;
                }
                
                uint8_t b1 = input ^ reg[5] ^ reg[3] ^ reg[2] ^ reg[0];
                uint8_t b2 = input ^ reg[3] ^ reg[2] ^ reg[1] ^ reg[0];
                
                expected_output[state][input] = (b1 << 1) | b2;
            }
        }
    }
    
    inline int getNextState(int current_state, int input_bit) const {
        return ((current_state << 1) | input_bit) & (NUM_STATES - 1);
    }
    
    inline int hammingDistance(uint8_t expected, uint8_t received) const {
        uint8_t diff = expected ^ received;
        return ((diff >> 1) & 1) + (diff & 1);
    }
    
public:
    ViterbiDecoder80211a() {
        computeExpectedOutput();
    }
    
    int decode(const int* encoded_bits, int n, int* output_bits) {
        if (n <= 0) return 0;
        
        // Матрица метрик
        int metrics[NUM_STATES];
        int next_metrics[NUM_STATES];
        
        // Матрица переходов для обратного прохода
        std::vector<std::vector<int>> transitions;
        transitions.resize(n + TAIL_BITS, std::vector<int>(NUM_STATES, -1));
        
        // Инициализация
        for (int i = 0; i < NUM_STATES; i++) {
            metrics[i] = INT_MAX;
        }
        metrics[0] = 0;
        
        // Шаг 1: Обработка реальных закодированных битов
        for (int step = 0; step < n; step++) {
            uint8_t received = (encoded_bits[2*step] << 1) | encoded_bits[2*step + 1];
            
            for (int i = 0; i < NUM_STATES; i++) {
                next_metrics[i] = INT_MAX;
            }
            
            for (int state = 0; state < NUM_STATES; state++) {
                if (metrics[state] == INT_MAX) continue;
                
                for (int input_bit = 0; input_bit < 2; input_bit++) {
                    int next_state = getNextState(state, input_bit);
                    uint8_t expected = expected_output[state][input_bit];
                    int distance = hammingDistance(expected, received);
                    int new_metric = metrics[state] + distance;
                    
                    if (new_metric < next_metrics[next_state]) {
                        next_metrics[next_state] = new_metric;
                        transitions[step][next_state] = (state << 1) | input_bit;
                    }
                }
            }
            
            std::copy(next_metrics, next_metrics + NUM_STATES, metrics);
        }
        
        // Шаг 2: Принудительные нулевые входные биты (хвост)
        for (int step = n; step < n + TAIL_BITS; step++) {
            for (int i = 0; i < NUM_STATES; i++) {
                next_metrics[i] = INT_MAX;
            }
            
            for (int state = 0; state < NUM_STATES; state++) {
                if (metrics[state] == INT_MAX) continue;
                
                int input_bit = 0;
                int next_state = getNextState(state, input_bit);
                uint8_t expected = expected_output[state][input_bit];
                
                int distance = 0;
                int new_metric = metrics[state] + distance;
                
                if (new_metric < next_metrics[next_state]) {
                    next_metrics[next_state] = new_metric;
                    transitions[step][next_state] = (state << 1) | input_bit;
                }
            }
            
            std::copy(next_metrics, next_metrics + NUM_STATES, metrics);
        }
        
        // Конечное состояние должно быть 0 (после хвостовых нулей)
        int final_state = 0;
        
        // Обратный проход
        std::vector<int> decoded_reversed;
        decoded_reversed.reserve(n + TAIL_BITS);
        
        int current_state = final_state;
        
        for (int step = n + TAIL_BITS - 1; step >= 0; step--) {
            int transition = transitions[step][current_state];
            if (transition == -1) {
                // Fallback: ищем любой переход
                for (int state = 0; state < NUM_STATES; state++) {
                    for (int bit = 0; bit < 2; bit++) {
                        if (getNextState(state, bit) == current_state) {
                            transition = (state << 1) | bit;
                            break;
                        }
                    }
                    if (transition != -1) break;
                }
            }
            
            int input_bit = transition & 1;
            decoded_reversed.push_back(input_bit);
            current_state = transition >> 1;
        }
        
        // Инвертируем и отбрасываем хвостовые биты
        int output_len = n;
        for (int i = 0; i < output_len; i++) {
            output_bits[i] = decoded_reversed[output_len + TAIL_BITS - 1 - i];
        }
        
        return output_len;
    }
};

static ViterbiDecoder80211a decoder;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("viterbi_decode:invalidInput",
                         "Need input argument.");
    }
    
    double* encoded_double = mxGetPr(prhs[0]);
    mwSize n_elements = mxGetNumberOfElements(prhs[0]);
    
    if (n_elements % 2 != 0) {
        mexErrMsgIdAndTxt("viterbi_decode:invalidInput",
                         "Input length must be even.");
    }
    
    int n = n_elements / 2;
    std::vector<int> encoded_bits(n_elements);
    for (mwSize i = 0; i < n_elements; i++) {
        encoded_bits[i] = (encoded_double[i] != 0) ? 1 : 0;
    }
    
    std::vector<int> output_bits(n);
    int output_len = decoder.decode(encoded_bits.data(), n, output_bits.data());
    
    plhs[0] = mxCreateDoubleMatrix(1, output_len, mxREAL);
    double* output = mxGetPr(plhs[0]);
    
    for (int i = 0; i < output_len; i++) {
        output[i] = output_bits[i];
    }
}

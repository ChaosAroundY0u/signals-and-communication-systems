#include <iostream>
#include <vector>
#include <cstdint>
#include "mex.h"

class ConvolutionalEncoder80211a {
private:
    // Сдвиговый регистр (6 бит памяти + текущий бит)
    uint8_t shift_register = 0;
    
    // Вычисление выходных битов
    void encodeBit(uint8_t input_bit, uint8_t& out_a, uint8_t& out_b) {
        // Сдвиг регистра и добавление нового бита
        shift_register = ((shift_register << 1) | input_bit) & 0x7F;
        
        // Выход A (полином 133) - исправленная версия
        out_a = (
            ((shift_register >> 6) & 1) ^ ((shift_register >> 4) & 1) ^
            ((shift_register >> 3) & 1) ^ ((shift_register >> 1) & 1) ^
            ((shift_register >> 0) & 1)
        );
        
        // Выход B (полином 171) - исправленная версия
        out_b = (
            ((shift_register >> 6) & 1) ^ ((shift_register >> 5) & 1) ^
            ((shift_register >> 4) & 1) ^ ((shift_register >> 3) & 1) ^
            ((shift_register >> 0) & 1)
        );
    }
    
public:
    void reset() {
        shift_register = 0;
    }
    
    // Кодирование последовательности битов (для MEX)
    // Вход: input_bits - массив из n битов (0/1)
    // Выход: output_bits - массив из 2*n битов (выделен caller'ом)
    void encode(const uint8_t* input_bits, int n, uint8_t* output_bits) {
        for (int i = 0; i < n; i++) {
            uint8_t a, b;
            encodeBit(input_bits[i], a, b);
            output_bits[2*i] = a;
            output_bits[2*i + 1] = b;
        }
    }
    
    // Кодирование с добавлением хвостовых нулей (6 битов для очистки регистра)
    // Вход: input_bits - массив из n битов (0/1)
    // Выход: output_bits - массив из 2*(n+6) битов (выделен caller'ом)
    void encodeWithTail(const uint8_t* input_bits, int n, uint8_t* output_bits) {
        // Кодируем исходные биты
        for (int i = 0; i < n; i++) {
            uint8_t a, b;
            encodeBit(input_bits[i], a, b);
            output_bits[2*i] = a;
            output_bits[2*i + 1] = b;
        }
        
        // Кодируем 6 хвостовых нулей для очистки регистра
        for (int i = 0; i < 6; i++) {
            uint8_t a, b;
            encodeBit(0, a, b);
            output_bits[2*(n + i)] = a;
            output_bits[2*(n + i) + 1] = b;
        }
    }
};


static ConvolutionalEncoder80211a encoder;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Проверка аргументов
    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("conv_encode:invalidInput",
                         "Need 1 or 2 input arguments.\n"
                         "Usage: encoded = conv_encode(bits, use_tail)");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("conv_encode:invalidOutput",
                         "Too many output arguments.");
    }
    
    // Получаем входной массив битов
    if (!mxIsDouble(prhs[0]) && !mxIsInt32(prhs[0]) && !mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("conv_encode:invalidInput",
                         "Input must be numeric or logical array.");
    }
    
    double* input_double = mxGetPr(prhs[0]);
    mwSize n_bits = mxGetNumberOfElements(prhs[0]);
    
    // Конвертируем в uint8_t массив
    std::vector<uint8_t> input_bits(n_bits);
    for (mwSize i = 0; i < n_bits; i++) {
        input_bits[i] = (input_double[i] != 0) ? 1 : 0;
    }
    
    // Параметр use_tail (опционально, по умолчанию 1 - добавлять хвостовые нули)
    int use_tail = 1;
    if (nrhs >= 2) {
        use_tail = (int)mxGetScalar(prhs[1]);
    }
    
    // Сбрасываем кодер перед кодированием
    encoder.reset();
    
    // Кодируем
    int output_len;
    std::vector<uint8_t> output_bits;
    
    if (use_tail) {
        // Добавляем 6 хвостовых нулей
        output_len = 2 * (n_bits + 6);
        output_bits.resize(output_len);
        encoder.encodeWithTail(input_bits.data(), n_bits, output_bits.data());
    } else {
        // Без хвостовых нулей
        output_len = 2 * n_bits;
        output_bits.resize(output_len);
        encoder.encode(input_bits.data(), n_bits, output_bits.data());
    }
    
    // Создаем выходной массив
    plhs[0] = mxCreateDoubleMatrix(1, output_len, mxREAL);
    double* output = mxGetPr(plhs[0]);
    
    for (int i = 0; i < output_len; i++) {
        output[i] = output_bits[i];
    }
}

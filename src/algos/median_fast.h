#pragma once

float median3x3(float* array);
float median24(float *array);
float median5x5(float* array);
float median7x7(float* array);
float median9x9(float* array);

inline __attribute__((always_inline)) float mymin(float a, float b) {
    return b < a ? b : a;
}

inline __attribute__((always_inline))float mymax(float a, float b) {
    return a < b ? b : a;
}

inline __attribute__((always_inline)) float median9f(float array0, float array1, float array2, float array3, float array4, float array5, float array6, float array7, float array8)
{
    float tmp = mymin(array1, array2);
    array2 = mymax(array1, array2);
    array1 = tmp;
    tmp = mymin(array4, array5);
    array5 = mymax(array4, array5);
    array4 = tmp;
    tmp = mymin(array7, array8);
    array8 = mymax(array7, array8);
    array7 = tmp;
    tmp = mymin(array0, array1);
    array1 = mymax(array0, array1);
    array0 = tmp;
    tmp = mymin(array3, array4);
    array4 = mymax(array3, array4);
    array3 = tmp;
    tmp = mymin(array6, array7);
    array7 = mymax(array6, array7);
    array6 = tmp;
    tmp = mymin(array1, array2);
    array2 = mymax(array1, array2);
    array1 = tmp;
    tmp = mymin(array4, array5);
    array5 = mymax(array4, array5);
    array4 = tmp;
    tmp = mymin(array7, array8);
    array8 = mymax(array7, array8);
    array3 = mymax(array0, array3);
    array5 = mymin(array5, array8);
    array7 = mymax(array4, tmp);
    tmp = mymin(array4, tmp);
    array6 = mymax(array3, array6);
    array4 = mymax(array1, tmp);
    array2 = mymin(array2, array5);
    array4 = mymin(array4, array7);
    tmp = mymin(array4, array2);
    array2 = mymax(array4, array2);
    array4 = mymax(array6, tmp);
    return mymin(array4, array2);
}

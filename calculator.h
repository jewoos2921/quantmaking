//
// Created by jewoo on 2022-10-05.
//

#pragma once

template<typename dataType>
class CCalculator {
public:
    CCalculator() = default;

    ~CCalculator() = default;

public:
    dataType addition(dataType a, dataType b);

    dataType subtraction(dataType a, dataType b);

    dataType multiplication(dataType a, dataType b);

    dataType division(dataType a, dataType b);
};

template<typename dataType>
dataType CCalculator<dataType>::addition(dataType a, dataType b) {
    return a + b;
}

template<typename dataType>
dataType CCalculator<dataType>::subtraction(dataType a, dataType b) {
    return a - b;
}

template<typename dataType>
dataType CCalculator<dataType>::multiplication(dataType a, dataType b) {
    return a * b;
}

template<typename dataType>
dataType CCalculator<dataType>::division(dataType a, dataType b) {
    return a / b;
}

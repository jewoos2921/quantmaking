#include <iostream>
#include "index.h"
#include "yield.h"
#include "option.h"

int main() {
    CIndex index;
    CYield yield;
    COption option;


    index.get_spot();
    index.get_vol();
    index.get_dividend();

    yield.get_riskfree();
    option.get_strike();
    option.get_maturity();
    option.european_calloption_price(index, yield);

    std::cout << "European Call option 가격: " << option.m_optionprice << std::endl;
    return 0;
}

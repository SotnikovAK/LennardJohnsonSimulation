#pragma warning(disable : 4996)
#pragma warning(disable : 4244)

#include "Libraries.h"
#include "Constants.h"

#include "Element.h"
#include "Simulation.h"

int main()
{

    int N;
    int BorderSize_;
    int flag_animation;


    std::cout << "Simulation Starts!" << std::endl;
    std::cout << "Enter the required number of elements" << std::endl;

    std::cin >> N;

    std::cout << "Enter size of space:" << std::endl;
    std::cin >> BorderSize_;

    std::cout << "Do I need to display the image? (1/0)" << std::endl;
    std::cout << "(at large values of n the animation lags a lot - for animation it is recommended up to 30)" << std::endl;
    std::cout << "(Animation is not for analysis - data is not displayed;)" << std::endl;
    std::cout << "Detailed settings in Constants.h" << std::endl;
    std::cin >> flag_animation;

    Simulation solution(BorderSize_, N, flag_animation);
    solution.WindowInteraction();

    return 0;
}
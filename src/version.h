
#ifndef SPINWALK_VERSION_H
#define SPINWALK_VERSION_H

#include <iostream>

#define SPINWALK_VERSION_MAJOR 1
#define SPINWALK_VERSION_MINOR 13
#define SPINWALK_VERSION_PATCH 10

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
inline void print_logo()
{ 
 std::cout << " \n"
" ____            _          __        __          _   _        \n"
"/ ___|   _ __   (_)  _ __   \\ \\      / /   __ _  | | | | __    \n"
"\\___ \\  | '_ \\  | | | '_ \\   \\ \\ /\\ / /   / _` | | | | |/ /    \n"
" ___) | | |_) | | | | | | |   \\ V  V /   | (_| | | | |   <     \n"
"|____/  | .__/  |_| |_| |_|    \\_/\\_/     \\__,_| |_| |_|\\_\\    \n"
"        |_|                                                    \n\n";

std::cout << "SpinWalk ver. " << SPINWALK_VERSION_MAJOR << "." << SPINWALK_VERSION_MINOR << "." << SPINWALK_VERSION_PATCH << std::endl;
}

#endif // SPINWALK_VERSION_H


#ifndef SPINWALK_VERSION_H
#define SPINWALK_VERSION_H

#include <iostream>
#include <string>

#define SPINWALK_VERSION_MAJOR 1
#define SPINWALK_VERSION_MINOR 15
#define SPINWALK_VERSION_PATCH 0

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
inline std::string get_verion()
{
  return std::to_string(SPINWALK_VERSION_MAJOR) + "." + std::to_string(SPINWALK_VERSION_MINOR) + "." + std::to_string(SPINWALK_VERSION_PATCH);
}

inline void print_logo()
{ 
 std::cout << " \n"
" ____            _          __        __          _   _        \n"
"/ ___|   _ __   (_)  _ __   \\ \\      / /   __ _  | | | | __    \n"
"\\___ \\  | '_ \\  | | | '_ \\   \\ \\ /\\ / /   / _` | | | | |/ /    \n"
" ___) | | |_) | | | | | | |   \\ V  V /   | (_| | | | |   <     \n"
"|____/  | .__/  |_| |_| |_|    \\_/\\_/     \\__,_| |_| |_|\\_\\    \n"
"        |_|                                                    \n\n";

std::cout << "SpinWalk ver. " << get_verion() << std::endl;
}



#endif // SPINWALK_VERSION_H

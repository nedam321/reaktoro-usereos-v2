#include <reaktoro-usereos.hpp>

#include <iostream>

int main() {
  std::cout << "Hello world!" << std::endl;
  int a[] = {0,0};
  int b[] = {2};
  int c[] = {0,0};
  GetDimensions(a,b,c);
  std::cout << "Hello world!" << std::endl;
  return 0;
}
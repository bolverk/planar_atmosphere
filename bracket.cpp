#include "bracket.hpp"

Bracket::Bracket(int low, int high):
  lh_(low, high) {}

bool Bracket::operator()(int num) const
{
  return num>=lh_.first && num<=lh_.second;
}

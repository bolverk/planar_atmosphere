#ifndef BRACKET_HPP
#define BRACKET_HPP 1

#include <utility>

using std::pair;

class Bracket
{
public:

  Bracket(int low, int high);

  bool operator()(int num) const;

private:

  const pair<int,int> lh_;
};

#endif // BRACKET_HPP

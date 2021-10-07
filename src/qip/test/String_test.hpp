#pragma once
#include "qip/String.hpp"
#include <cassert>

void String_test() {
  assert(qip::wildcard_compare("helloben", "helloben") == true);
  assert(qip::wildcard_compare("helloben", "*ben") == true);
  assert(qip::wildcard_compare("helloben", "hello*") == true);
  assert(qip::wildcard_compare("helloben", "hello*ben") == true);
  assert(qip::wildcard_compare("helloben", "hell*ben") == true);
  assert(qip::wildcard_compare("helloben", "he*en") == true);
  assert(qip::wildcard_compare("helloben", "*") == true);
  assert(qip::wildcard_compare("helloben", "*helloben") == true);
  assert(qip::wildcard_compare("helloben", "helloben*") == true);
  assert(qip::wildcard_compare("", "*") == true);
  assert(qip::wildcard_compare("", "") == true);

  assert(qip::wildcard_compare("helloben", "xhelloben") == false);
  assert(qip::wildcard_compare("helloben", "x*ben") == false);
  assert(qip::wildcard_compare("helloben", "*benx") == false);
  assert(qip::wildcard_compare("helloben", "xhello*") == false);
  assert(qip::wildcard_compare("helloben", "hello*x") == false);
  assert(qip::wildcard_compare("helloben", "*x") == false);
  assert(qip::wildcard_compare("helloben", "x*") == false);
  assert(qip::wildcard_compare("helloben", "**") == false);
  assert(qip::wildcard_compare("helloben", "") == false);
  assert(qip::wildcard_compare("", "helloben") == false);
}

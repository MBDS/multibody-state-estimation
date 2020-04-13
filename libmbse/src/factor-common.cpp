/*+-------------------------------------------------------------------------+
  |    FactorGraph Control (fgcontrol)  C++ library                         |
  |                                                                         |
  | Copyright (C) 2019-2020 University of Almeria                           |
  | See README for list of authors and papers                               |
  | Distributed under GNU General Public License version 3                  |
  |   See <http://www.gnu.org/licenses/>                                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factor-common.h>
#include <cmath>
#include <iostream>

using namespace mbse;

void state_t::print(const std::string& prefix) const
{
	std::cout << prefix << this->transpose();
}

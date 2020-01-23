#include "extrap.hh"
#include "level++.hh"

template class traverse<tra_positive,ex_ref_map>;
//template class traverse<tra_negative,ex_ref_map>;
template class traverse<tra_positive,ex_staggered>;
//template class traverse<tra_negative,ex_staggered>;
template class traverse<tra_positive,ex_regular>;
//template class traverse<tra_negative,ex_regular>;

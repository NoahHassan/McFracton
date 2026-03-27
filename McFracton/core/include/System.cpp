#include "System.h"

#include <assert.h>

System::~System() {}

double System::getSite(int index) const
{
    assert(0 <= index);
    assert(index < site_fields.size());
    return site_fields[index];
}

double System::getPlaq(int index) const
{
    assert(0 <= index);
    assert(index < plaq_fields.size());
    return plaq_fields[index];
}

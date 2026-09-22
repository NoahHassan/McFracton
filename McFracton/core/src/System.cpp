#include "System.h"

#include <assert.h>

System::~System() {}

double System::variable(int index) const
{
    assert(0 <= index);
    assert(index < fields.size());
    return fields[index];
}

void System::applyUpdate(int index, double delta)
{
    assert(0 <= index);
    assert(index < fields.size());
    fields[index] += delta;
}

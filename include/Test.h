#ifndef IMPBFF_TEST_H
#define IMPBFF_TEST_H

#include <IMP/bff/bff_config.h>
#include <IMP/Object.h>
#include <IMP/object_macros.h>
#include <IMP/warning_macros.h>

#include <iostream>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT BFFTest : public IMP::Object {

private:

public:

    void show();

    BFFTest() : IMP::Object("BFFTest%1%") {}

    IMP_OBJECT_METHODS(BFFTest);
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_TEST_H
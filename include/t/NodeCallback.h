#ifndef chinet_NODECALLBACK_H
#define chinet_NODECALLBACK_H

#include <IMP/bff/bff_config.h>
#include <IMP/Object.h>
#include <IMP/object_macros.h>
#include <IMP/warning_macros.h>

#include <functional>
#include <algorithm> /* std::min std::max */
#include <rttr/registration>

#include "Port.h"
#include "CNode.h"
#include "../Functions.h"

class Port;


IMPBFF_BEGIN_NAMESPACE

class NodeCallback{

public:
    
    virtual void run(
            std::map<std::string, std::shared_ptr<Port>>,
            std::map<std::string, std::shared_ptr<Port>>
            ){
        std::cout << "This print by NodeCallback class from C++" << std::endl;
    }

    NodeCallback() = default;
    virtual ~NodeCallback() {};
};

IMPBFF_END_NAMESPACE

#endif //chinet_NODECALLBACK_H

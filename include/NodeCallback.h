#ifndef IMPBFF_NODECALLBACK_H
#define IMPBFF_NODECALLBACK_H

#include <IMP/bff/bff_config.h>

#include <map>
#include <string>
#include <functional>
#include <algorithm> /* std::min std::max */
#include <rttr/registration>


IMPBFF_BEGIN_NAMESPACE
class Port;

#include <IMP/bff/Port.h>
#include <IMP/bff/Node.h>
#include <IMP/bff/internal/Functions.h>


class IMPBFFEXPORT NodeCallback{

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

#endif //IMPBFF_NODECALLBACK_H

#ifndef IMPBFF_NODECALLBACK_H
#define IMPBFF_NODECALLBACK_H

#include <IMP/bff/bff_config.h>

#include <map>
#include <string>
#include <functional>
#include <algorithm> /* std::min std::max */


IMPBFF_BEGIN_NAMESPACE
class CnPort;

#include <IMP/bff/CnPort.h>
#include <IMP/bff/CnNode.h>
#include <IMP/bff/internal/Functions.h>



class IMPBFFEXPORT CnNodeCallback{

public:
    
    virtual void run(
            std::map<std::string, std::shared_ptr<CnPort>>,
            std::map<std::string, std::shared_ptr<CnPort>>
            ){
        std::cout << "This print by NodeCallback class from C++" << std::endl;
    }

    CnNodeCallback() = default;
    virtual ~CnNodeCallback() {};
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_NODECALLBACK_H

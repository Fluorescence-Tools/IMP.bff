#include <IMP/bff/NodeCallback.h>

IMPBFF_BEGIN_NAMESPACE

//using namespace rttr;

template <typename T>
inline void mul(T* tmp, size_t &n_elements, const std::map<std::string, std::shared_ptr<Port>> &inputs){
    std::fill(tmp, tmp + n_elements, 1.0);
    for(auto &o : inputs){
        T* va; int vn;
        o.second->get_value<T>(&va, &vn);
        for(size_t i=0; i < n_elements; i++){
            tmp[i] *= va[i];
        }
    }
}

template <typename T>
inline void add(
        T* tmp,
        size_t &n_elements,
        const std::map<std::string, std::shared_ptr<Port>> &inputs
)
{
    std::fill(tmp, tmp + n_elements, 1.0);
    for(auto &o : inputs){
        T* va; int vn;
        o.second->get_value<T>(&va, &vn);
        for(size_t i=0; i < n_elements; i++){
            tmp[i] += va[i];
        }
    }
}


template <typename T>
void combine(
        std::map<std::string, std::shared_ptr<Port>> &inputs,
        std::map<std::string, std::shared_ptr<Port>> &outputs,
        int operation
        ){
#if IMPBFF_VERBOSE
    std::clog << "-- Combining values of input ports."  << std::endl;
#endif
    size_t n_elements = UINT_MAX;
#if IMPBFF_VERBOSE
    std::clog << "-- Determining input vector with smallest length." << std::endl;
#endif
    for(auto &o : inputs){
        n_elements = std::min(n_elements, o.second->size());
    }
#if IMPBFF_VERBOSE
    std::clog << "-- The smallest input vector has a length of: " << n_elements << std::endl;
#endif
    auto tmp = (T*) malloc(n_elements * sizeof(T));
    switch(operation){
        case BFF_PORT_ADD:
#if IMPBFF_VERBOSE
            std::clog << "-- Adding input ports" << std::endl;
#endif
            add(tmp, n_elements, inputs);
            break;
        case BFF_PORT_MUL:
#if IMPBFF_VERBOSE
            std::clog << "-- Multiplying input ports" << std::endl;
#endif
            mul(tmp, n_elements, inputs);
            break;
        default:
            break;
    }
    if(!outputs.empty()){
        if (outputs.find("outA") == outputs.end() ) {
#if IMPBFF_VERBOSE
            std::clog << "ERROR: Node does not define output port with the name 'outA' " << std::endl;
#endif
            outputs.begin()->second->set_value(tmp, n_elements);
        } else {
#if IMPBFF_VERBOSE
            std::clog << "Setting value to output " << std::endl;
#endif
            outputs["outA"]->set_value(tmp, n_elements);
        }
    } else{
        std::cerr << "ERROR: There are no output ports." << std::endl;
    }
    free(tmp);
}


template <typename T>
void addition(
        std::map<std::string, std::shared_ptr<Port>> &inputs,
        std::map<std::string, std::shared_ptr<Port>> &outputs
)
{
#if IMPBFF_VERBOSE
    std::clog << "addition"  << std::endl;
#endif
    combine<T>(inputs, outputs, BFF_PORT_ADD);
}

template <typename T>
void multiply(
        std::map<std::string, std::shared_ptr<Port>> &inputs,
        std::map<std::string, std::shared_ptr<Port>> &outputs
)
{
    combine<T>(inputs, outputs, BFF_PORT_MUL);
}

void nothing(
        std::map<std::string, std::shared_ptr<Port>> &inputs,
        std::map<std::string, std::shared_ptr<Port>> &outputs){
}

void passthrough(
        std::map<std::string, std::shared_ptr<Port>> &inputs,
        std::map<std::string, std::shared_ptr<Port>> &outputs
        ){
    for(auto it_in = inputs.cbegin(), end_in = inputs.cend(),
            it_out = outputs.cbegin(), end_out = outputs.cend();
            it_in != end_in || it_out != end_out;)
    {
        double* va; int vn;
        it_in->second->get_value(&va, &vn);
        auto v = std::vector<double>();
        v.assign(va, va + vn);
        if(it_in != end_in) {
            ++it_in;
        } else{
            break;
        }
        if(it_out != end_out) {
            it_out->second->set_value(v.data(), v.size());
            ++it_out;
        } else{
            break;
        }
    }
}


RTTR_REGISTRATION{
    using namespace rttr;
    registration::class_<NodeCallback>("NodeCallback").constructor<>().method("run", &NodeCallback::run);
    registration::method("multiply_double", &multiply<double>);
    registration::method("multiply_int", &multiply<long>);
    registration::method("addition_double", &addition<double>);
    registration::method("addition_int", &addition<long>);
    registration::method("nothing", &nothing);
    registration::method("passthrough", &passthrough);
}

IMPBFF_END_NAMESPACE

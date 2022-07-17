#include <IMP/bff/bff_config.h>
#include <IMP/bff/Port.h>

IMPBFF_BEGIN_NAMESPACE


std::vector<unsigned char> Port::get_bytes(){
    std::vector<unsigned char> output;
    size_t n_output = n_buffer_elements_ * buffer_element_size_;
    unsigned char *charBuffer = (unsigned char *) buffer_;
    output.insert(output.end(), charBuffer[0], charBuffer[n_output]);
    return output;
}

void Port::set_bytes(unsigned char *input, int n_input){
    buffer_ = realloc(buffer_, n_input);
    if(buffer_ != nullptr){
        memcpy(buffer_, input, n_input);
    }
}

// Methods
// --------------------------------------------------------------------
bool Port::read_from_db(const std::string &oid_string)
{
    bool re = MongoObject::read_from_db(oid_string);
    auto v = MongoObject::get_array<double>("value");
    n_buffer_elements_ = v.size();
    buffer_ = std::realloc(buffer_, sizeof(double) * n_buffer_elements_);
    if(buffer_ != nullptr)
    {
        n_buffer_elements_ = v.size();
        if(value_type == 0){
            memcpy(buffer_, v.data(), n_buffer_elements_ * sizeof(double));
        } else{
            memcpy(buffer_, v.data(), n_buffer_elements_ * sizeof(int));
        }
    }
    return re;
}

bool Port::write_to_db()
{
    bson_t doc = get_bson();
    return MongoObject::write_to_db(doc, 0);
}

// Setter
//--------------------------------------------------------------------

bson_t Port::get_bson()
{
    bson_t dst = get_bson_excluding("value", "bounds", NULL);
    if((value_type == BFF_PORT_VECTOR_FLOAT) || (value_type == BFF_PORT_FLOAT)){
        long* va; int nv;
        get_own_value(&va, &nv);
        auto v = std::vector<double>();
        v.assign(va, va + nv);
        append_number_array(&dst, "value", v);
    } else{
        double* va; int nv;
        get_own_value(&va, &nv);
        auto v = std::vector<double>();
        v.assign(va, va + nv);
        append_number_array(&dst, "value", v);
    }
    append_number_array(&dst, "bounds", bounds_);
    return dst;
}

bool Port::bound_is_valid()
{
    if (bounds_.size() == 2) {
        if (bounds_[0] != bounds_[1]) {
            return true;
        }
    }
    return false;
}


void Port::set_bounds(std::vector<double> b)
{
    if (b.size() >= 2) {
        bounds_.clear();
        bounds_.push_back(std::min(b[0], b[1]));
        bounds_.push_back(std::max(b[0], b[1]));
    }
}

std::vector<double> Port::get_bounds()
{
    return bounds_;
}




// Getter
//--------------------------------------------------------------------


bool Port::is_fixed()
{
    return get_<bool>("is_fixed");
}

void Port::set_fixed(bool fixed)
{
    set_("is_fixed", fixed);
}

bool Port::is_output()
{
    return get_<bool>("is_output");
}

void Port::set_port_type(bool is_output)
{
    set_("is_output", is_output);
}

bool Port::is_reactive()
{
    return get_<bool>("is_reactive");
}

void Port::set_reactive(bool reactive)
{
    set_("is_reactive", reactive);
}

bool Port::is_bounded()
{
    return get_<bool>("is_bounded");
}

void Port::set_bounded(bool bounded) {
    set_("is_bounded", bounded);
}

std::shared_ptr<Node> Port::get_node() {
    return node_;
}

void Port::set_node(std::shared_ptr<Node> node_ptr) {
    node_ = node_ptr;
}

bool Port::remove_links_to_port() {
    if (link_ != nullptr) {
        // remove pointer to this port in the port to which this is linked
        auto it = std::find(linked_to_.begin(), linked_to_.end(), this);
        if (it != linked_to_.end()) {
            linked_to_.erase(it);
            return true;
        }
    }
    return false;
}

void Port::update_attached_node() {
#if IMPBFF_VERBOSE
    std::clog << "-- Port is reactive: " << is_reactive() << std::endl;
    std::clog << "-- Port is output: " << is_output() << std::endl;
    std::clog << "-- Updating the node: " << node_->get_name() << std::endl;
#endif
    node_->set_valid(false);
    if (is_reactive() && !is_output()) {
        node_->evaluate();
    }
#if IMPBFF_VERBOSE
    std::clog << "-- Node is valid: " << node_->is_valid() << std::endl;
#endif
}

IMPBFF_END_NAMESPACE

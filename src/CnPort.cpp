#include <IMP/bff/CnPort.h>

IMPBFF_BEGIN_NAMESPACE


CnPort::CnPort(
        bool fixed,
        bool is_output,
        bool is_reactive,
        bool is_bounded,
        double lb,
        double ub,
        int value_type,
        std::string name
) : CnMongoObject(name) {

    append_string(&document, "type", "port");
    bson_append_bool(&document, "is_output", 9, false);
    bson_append_bool(&document, "is_fixed", 8, false);
    bson_append_bool(&document, "is_reactive", 11, false);
    bson_append_oid(&document, "link", 4, &oid_document);
    bson_append_bool(&document, "is_bounded", 10, false);

    set_fixed(fixed);
    set_port_type(is_output);
    set_reactive(is_reactive);
    set_bounded(is_bounded);

    if (is_bounded) {
        bounds_.push_back(lb);
        bounds_.push_back(ub);
    }

    CnPort::value_type = value_type;
}

void CnPort::get_bytes(unsigned char **output, int *n_output, bool copy){
    *n_output = cc_.size();
    if(copy){
        auto out = malloc(cc_.size());
        memcpy(out, cc_.get_bytes().data(), cc_.size());
        *output = static_cast<unsigned char*>(out);
    } else{
        *output = cc_.get_bytes().data();
    }
}

void CnPort::set_bytes(unsigned char *input, int n_input){
    cc_.set_bytes(input, n_input);
}

// Methods
// --------------------------------------------------------------------
bool CnPort::read_from_db(const std::string &oid_string)
{
    bool re = CnMongoObject::read_from_db(oid_string);
    auto v = CnMongoObject::get_array<unsigned char>("value");
    cc_.set_bytes(v);
    return re;
}

bool CnPort::write_to_db()
{
    bson_t doc = get_bson();
    return CnMongoObject::write_to_db(doc, 0);
}

// Setter
//--------------------------------------------------------------------

bson_t CnPort::get_bson()
{
    bson_t dst = get_bson_excluding("value", "bounds", NULL);
    if(value_type == 0){
        long* va; int nv;
        get_own_value(&va, &nv);
        auto v = std::vector<long>();
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

bool CnPort::bound_is_valid()
{
    if (bounds_.size() == 2) {
        if (bounds_[0] != bounds_[1]) {
            return true;
        }
    }
    return false;
}


void CnPort::set_bounds(double *input, int n_input)
{
    if (n_input >= 2) {
        bounds_.clear();
        double lower = std::min(input[0], input[1]);
        double upper = std::max(input[0], input[1]);
        bounds_.push_back(lower);
        bounds_.push_back(upper);
    }
}

void CnPort::get_bounds(double **output, int *n_output)
{
    *output = bounds_.data();
    *n_output = bounds_.size();
}

// Getter
//--------------------------------------------------------------------

bool CnPort::is_fixed()
{
    return get_<bool>("is_fixed");
}

void CnPort::set_fixed(bool fixed)
{
    set_("is_fixed", fixed);
}

bool CnPort::is_output()
{
    return get_<bool>("is_output");
}

void CnPort::set_port_type(bool is_output)
{
    set_("is_output", is_output);
}

bool CnPort::is_reactive()
{
    return get_<bool>("is_reactive");
}

void CnPort::set_reactive(bool reactive)
{
    set_("is_reactive", reactive);
}

bool CnPort::is_bounded()
{
    return get_<bool>("is_bounded");
}

void CnPort::set_bounded(bool bounded)
{
    set_("is_bounded", bounded);
}

CnNode* CnPort::get_node()
{
    return node_;
}

void CnPort::set_node(CnNode *node_ptr)
{
    node_ = node_ptr;
}

bool CnPort::is_float() {
    return ((get_value_type() == BFF_PORT_VECTOR_FLOAT) ||
            (get_value_type() == BFF_PORT_FLOAT));
}

int CnPort::get_value_type() {
    if (!is_linked()) {
        return value_type;
    } else {
        return get_link()->value_type;
    }
}

void CnPort::set_value_type(int v) {
    value_type = v;
    if (is_linked()) {
        get_link()->set_value_type(v);
    }
}

void CnPort::set_link(std::shared_ptr<CnPort> v) {
#if IMPBFF_VERBOSE
    std::clog << "SET_LINK" << std::endl;
std::clog << "-- Link to Port: " << v->get_name() << std::endl;
std::clog << "-- Link value type: " << v->value_type << std::endl;
#endif
    if (v != nullptr) {
        if(v.get() == this){
#if IMPBFF_VERBOSE
            std::clog << "WARNING: cannot link to self." << std::endl;
#endif
        } else{
            set_oid("link", v->get_bson_oid());
            link_ = v;
            v->linked_to_.push_back(this);
        }
    } else {
        unlink();
    }
    if(node_!=nullptr) update_attached_node();
}

bool CnPort::unlink() {
    set_oid("link", get_bson_oid());
    bool re = true;
    re &= remove_links_to_port();
    link_ = nullptr;
    return re;
}

bool CnPort::is_linked() {
    return (link_ != nullptr);
}

std::vector<CnPort *> CnPort::get_linked_ports() {
    return linked_to_;
}

std::shared_ptr<CnPort> CnPort::get_link() {
    return link_;
}

void CnPort::update_attached_node() {
#if IMPBFF_VERBOSE
    std::clog << "-- CnPort is reactive: " << is_reactive() << std::endl;
    std::clog << "-- CnPort is output: " << is_output() << std::endl;
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

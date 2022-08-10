#include <IMP/bff/CnNode.h>

IMPBFF_BEGIN_NAMESPACE

CnNode::CnNode(
        std::string name,
        const std::map<std::string, std::shared_ptr<CnPort>>& ports,
        std::shared_ptr<CnNodeCallback>
) : CnMongoObject(name)
{
    append_string(&document, "type", "CnNode");
    set_ports(ports);
    this->callback_class = callback_class;
}


// Destructor
//--------------------------------------------------------------------
CnNode::~CnNode() = default;


// Methods
//--------------------------------------------------------------------

bool CnNode::read_from_db(const std::string &oid_string){
#if IMPBFF_VERBOSE
    std::clog << "TODO!! READING CnNode FROM DB" << std::endl;
    std::clog << "Requested OID:" << oid_string << std::endl;
#endif
    bool return_value = true;
    return_value &= CnMongoObject::read_from_db(oid_string);
    return_value &= create_and_connect_objects_from_oid_doc(
            &document, "ports", &ports);
#if IMPBFF_VERBOSE
    std::clog << "callback-restore: " << get_string_by_key(&document, "callback") << std::endl;
#endif
    return return_value;
}

bool CnNode::write_to_db() {
    bool re = CnMongoObject::write_to_db();
    for(auto &o : ports){
        if(!o.second->is_connected_to_db()){
            re &= connect_object_to_db(o.second);
        }
        o.second->write_to_db();
    }

    return re;
}

// Getter
//--------------------------------------------------------------------


//std::string CnNode::get_name(){
//    std::string r;
//    r.append(IMP::bff::CnMongoObject::get_name());
//    r.append("(");
//    for(auto const &n : get_input_ports()){
//        r.append(n.first);
//        r.append(",");
//    }
//    r.append(")");
//    r.append("->");
//    r.append("(");
//    for(auto const &n : get_output_ports()){
//        r.append(n.first);
//        r.append(",");
//    }
//    r.append(")");
//    return r;
//}


std::map<std::string, std::shared_ptr<CnPort>> CnNode::get_ports(){
    return ports;
}

void CnNode::set_ports(const std::map<std::string, std::shared_ptr<CnPort>>& ports){
    for(auto &o: ports){
        o.second->set_name(o.first);
        add_port(o.first, o.second, o.second->is_output(), false);
    }
    fill_input_output_port_lookups();
}

CnPort* CnNode::get_port(const std::string &port_name){
    return ports[port_name].get();
}

CnPort* CnNode::get_input_port(const std::string &port_name){
    return in_[port_name].get();
}

CnPort* CnNode::get_output_port(const std::string &port_name){
    return out_[port_name].get();
}

std::map<std::string, std::shared_ptr<CnPort>> CnNode::get_input_ports(){
    return in_;
}

std::map<std::string, std::shared_ptr<CnPort>> CnNode::get_output_ports(){
    return out_;
}


void CnNode::set_callback(std::shared_ptr<CnNodeCallback> cb){
    callback_class = cb;
}

void CnNode::add_port(
        const std::string &key,
        std::shared_ptr<CnPort> port, bool is_output, bool fill_in_out) {
#if IMPBFF_VERBOSE
    std::clog << "ADDING PORT TO CnNode" << std::endl;
    std::clog << "-- Name of CnNode: " << get_name() << std::endl;
    std::clog << "-- Key of port: " << key << std::endl;
    std::clog << "-- Port is_output: " << is_output << std::endl;
    std::clog << "-- Fill value of output: " << fill_in_out << std::endl;
#endif
    port->set_port_type(is_output);
    // auto n = std::dynamic_pointer_cast<CnNode>(MongoObject::shared_from_this());
    port->set_node(this);
    if (ports.find(key) == ports.end() ) {
#if IMPBFF_VERBOSE
        std::clog << "-- The key of the port was not found." << std::endl;
        std::clog << "-- Port " << key << " was created in CnNode. " << std::endl;
#endif
        ports[key] = port;
    } else {
        auto p = ports[key];
        if(port != p){
#if IMPBFF_VERBOSE
            std::clog << "WARNING: Overwriting the port that was originally associated to the key " << key << "." << std::endl;
#endif
            ports[key] = port;
        } else{
            std::cerr << "WARNING: Port is already part of the CnNode." << std::endl;
            std::cerr << "-- Assigning Port to the key: " << key << "." << std::endl;
        }
    }
    if(fill_in_out){
        fill_input_output_port_lookups();
    }
}

void CnNode::add_input_port(const std::string &key, std::shared_ptr<CnPort> port) {
    add_port(key, port, false, true);
}

void CnNode::add_output_port(const std::string &key, std::shared_ptr<CnPort> port) {
    add_port(key, port, true, true);
}

bson_t CnNode::get_bson(){
    bson_t dst = CnMongoObject::get_bson_excluding(
            "input_ports",
            "output_ports",
             "callback",
             "callback_type",
             NULL
    );
    create_oid_dict_in_doc<CnPort>(&dst, "ports", ports);
    return dst;
}

void CnNode::evaluate(){
#if IMPBFF_VERBOSE
    std::clog << "CnNode EVALUATE" << std::endl;
    std::clog << "-- CnNode name: " << get_name() << std::endl;
#endif
    callback_class->run(in_, out_);
#if IMPBFF_VERBOSE
    std::clog << "-- Setting CnNodes associated to output ports to invalid."  << std::endl;
#endif
    for(auto &o : get_output_ports()){
        auto n = o.second->get_node();
        if(n != nullptr){
#if IMPBFF_VERBOSE
            std::clog << "-- CnNode " << n->get_name() << " of port " << o.second->get_name() << " set to invalid." << std::endl;
#endif
            n->set_valid(false);
        }
    }
    node_valid_ = true;
}

void CnNode::fill_input_output_port_lookups(){
    out_.clear();
    in_.clear();
    for(auto &o: ports){
        if(o.second->is_output()){
            out_[o.first] = o.second;
        } else{
            in_[o.first] = o.second;
        }
    }
}

bool CnNode::inputs_valid(){
    for(const auto &i : in_){
        auto input_port = i.second;
        if(input_port->is_linked()){
            auto output_port = input_port->get_link();
            auto output_node = output_port->get_node();
            if(output_node == this) return true;
            else if(!output_node->is_valid()) return false;
        }
    }
    return true;
}

void CnNode::set_valid(bool is_valid){
    node_valid_ = is_valid;
    for(auto &v : out_)
    {
        auto output_port = v.second;
        //v.second->set_invalid();
    }
}

bool CnNode::is_valid(){
    if(get_input_ports().empty()) return true;
    else if(!inputs_valid()) return false;
    else return node_valid_;
}

IMPBFF_END_NAMESPACE

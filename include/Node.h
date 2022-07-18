#ifndef IMPBFF_NODE_H
#define IMPBFF_NODE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>

#include <IMP/bff/MongoObject.h>
#include <IMP/bff/Port.h>
#include <IMP/bff/NodeCallback.h>

IMPBFF_BEGIN_NAMESPACE

class Port;
class NodeCallback;
class IMPBFFEXPORT Node : public MongoObject{

private:
    friend Port;
    bool node_valid_ = false;
    std::map<std::string, std::shared_ptr<Port>> ports;

protected:
    std::map<std::string, std::shared_ptr<Port>> in_;
    std::map<std::string, std::shared_ptr<Port>> out_;

public:

    void fill_input_output_port_lookups();

    /// A pointer to a function that operates on an input Port instance (first argument)
    /// and writes to an output Port instance (second argument)
    std::shared_ptr<NodeCallback> callback_class;

    // Constructor & Destructor
    //--------------------------------------------------------------------
    Node(std::string name="",
         const std::map<std::string, std::shared_ptr<Port>>& ports =
                 std::map<std::string, std::shared_ptr<Port>>(),
         std::shared_ptr<NodeCallback> callback_class = nullptr
    );
    ~Node();

    // Methods
    //--------------------------------------------------------------------
    bool read_from_db(const std::string &oid_string) final;
    void evaluate();
    bool is_valid();
    bool inputs_valid();
    bool write_to_db() final;

    // Getter
    //--------------------------------------------------------------------
    bson_t get_bson();
    std::string get_name();
    std::map<std::string, std::shared_ptr<Port>> get_input_ports();
    std::map<std::string, std::shared_ptr<Port>> get_output_ports();
    std::map<std::string, std::shared_ptr<Port>> get_ports();
    void set_ports(const std::map<std::string, std::shared_ptr<Port>>& ports);
    Port* get_port(const std::string &port_name);
    void add_port(const std::string &key, std::shared_ptr<Port>,
            bool is_source, bool fill_in_out=true);
    void add_input_port(const std::string &key, std::shared_ptr<Port>);
    void add_output_port(const std::string &key, std::shared_ptr<Port>);
    Port* get_input_port(const std::string &port_name);
    Port* get_output_port(const std::string &port_name);

    // Setter
    //-----------------------------------------------------------
    void set_callback(std::shared_ptr<NodeCallback> cb);
    void set_valid(bool is_valid=false);
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_NODE_H

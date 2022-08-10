#ifndef IMPBFF_NODE_H
#define IMPBFF_NODE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>

#include <IMP/bff/CnMongoObject.h>
#include <IMP/bff/CnPort.h>
#include <IMP/bff/CnNodeCallback.h>

IMPBFF_BEGIN_NAMESPACE

class CnPort;
class CnNodeCallback;
class IMPBFFEXPORT CnNode : public CnMongoObject{

private:
    friend CnPort;
    bool node_valid_ = false;
    std::map<std::string, std::shared_ptr<CnPort>> ports;

protected:
    std::map<std::string, std::shared_ptr<CnPort>> in_;
    std::map<std::string, std::shared_ptr<CnPort>> out_;

public:

    void fill_input_output_port_lookups();

    /// A pointer to a function that operates on an input CnPort instance (first argument)
    /// and writes to an output CnPort instance (second argument)
    std::shared_ptr<CnNodeCallback> callback_class;

    // Constructor & Destructor
    //--------------------------------------------------------------------
    CnNode(std::string name="",
         const std::map<std::string, std::shared_ptr<CnPort>>& ports =
                 std::map<std::string, std::shared_ptr<CnPort>>(),
         std::shared_ptr<CnNodeCallback> callback_class = nullptr
    );
    ~CnNode();

    // Methods
    //--------------------------------------------------------------------
    bool read_from_db(const std::string &oid_string) final;
    void evaluate();
    bool is_valid();
    bool inputs_valid();
    bool write_to_db() final;

    // Getter / Setter
    //--------------------------------------------------------------------
    bson_t get_bson();
//    std::string get_name();
    std::map<std::string, std::shared_ptr<CnPort>> get_input_ports();
    std::map<std::string, std::shared_ptr<CnPort>> get_output_ports();
    std::map<std::string, std::shared_ptr<CnPort>> get_ports();
    void set_ports(const std::map<std::string, std::shared_ptr<CnPort>>& ports);
    void add_port(const std::string &key, std::shared_ptr<CnPort>,
            bool is_source, bool fill_in_out=true);
    void add_input_port(const std::string &key, std::shared_ptr<CnPort>);
    void add_output_port(const std::string &key, std::shared_ptr<CnPort>);
    CnPort* get_port(const std::string &port_name);
    CnPort* get_input_port(const std::string &port_name);
    CnPort* get_output_port(const std::string &port_name);

    // Setter
    //-----------------------------------------------------------
    void set_callback(std::shared_ptr<CnNodeCallback> cb);
    void set_valid(bool is_valid=false);
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_NODE_H

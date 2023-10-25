#ifndef IMPBFF_CNPORT_H
#define IMPBFF_CNPORT_H

#include <IMP/bff/bff_config.h>
#include <IMP/object_macros.h>

#include <cstdint>
#include <memory>
#include <algorithm> /* std::max, std::min */
#include <cmath>
#include <cstdlib>
#include <vector>
#include <map>
#include <utility> /* pair */
#include <string>

#include <IMP/bff/CnMongoObject.h>
#include <IMP/bff/CnNode.h>
#include <IMP/bff/internal/CnPortCache.h>
#include <IMP/bff/internal/Functions.h>

IMPBFF_BEGIN_NAMESPACE
class CnNode;

/// Specifies the type of the Port
typedef enum{
    BFF_PORT_VECTOR_INT,
    BFF_PORT_VECTOR_FLOAT,
    BFF_PORT_INT,
    BFF_PORT_FLOAT,
    BFF_PORT_BINARY
} CnPortType;


typedef enum{
    BFF_PORT_ADD,
    BFF_PORT_MUL
} PortOperation;



class IMPBFFEXPORT CnPort : public CnMongoObject {

private:

    std::map<int, size_t> PORT_VALUE_BYTE_SIZE = {
            {BFF_PORT_VECTOR_INT, sizeof(long)},
            {BFF_PORT_VECTOR_FLOAT, sizeof(double)},
            {BFF_PORT_INT, sizeof(long)},
            {BFF_PORT_FLOAT, sizeof(double)},
            {BFF_PORT_BINARY, 1}
    };

    CnPortCache cc_;
    std::vector<double> bounds_{};
    CnNode* node_ = nullptr;

    /*!
     * @brief Pointer to another Port (default value nullptr).
     * If not nullptr @class Port::get_value returns value of link
     */
    std::shared_ptr<CnPort> link_ = nullptr;

    /*!
     * @brief This attribute stores the Ports that are dependent on the value
     * of this Port object. If this Port object is reactive a change of the
     * value of this Port object is propagated to the dependent Ports.
     */
    std::vector<CnPort *> linked_to_;

    bool remove_links_to_port() {
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

    template<typename T>
    void set_value_of_dependents(T *input, int n_input) {
        for (auto &v : linked_to_) {
            v->set_value(input, n_input);
        }
    }

    /// Specifies the type of the Port
    int value_type = BFF_PORT_FLOAT;

    /// size (in bytes) of a port value
    int value_size = PORT_VALUE_BYTE_SIZE[BFF_PORT_FLOAT];


public:

    size_t size() {
        return cc_.size() / value_size;
    }

    // Constructor & Destructor
    //--------------------------------------------------------------------
    ~CnPort() {
        remove_links_to_port();
    };

    CnPort(
        bool fixed = false,
        bool is_output = false,
        bool is_reactive = false,
        bool is_bounded = false,
        double lb = 0,
        double ub = 0,
        int value_type = 1,
        std::string name = ""
    );

    void set_node(CnNode *node_ptr);

    CnNode* get_node();

    // Getter & Setter
    //--------------------------------------------------------------------
    int get_value_type();

    void set_value_type(int v);

    template<typename T>
    void update_value_type(int n_input){
        if (std::is_same<T, double>::value) {
            value_type = BFF_PORT_VECTOR_FLOAT;
        } else {
            value_type = BFF_PORT_VECTOR_INT;
        }
        if(n_input == 1)
            value_type += 2; // use scalar type for single numbers
    }

    template<typename T>
    void set_value(T *input, int n_input) {
#if IMPBFF_VERBOSE
        std::clog << "SET PORT VALUE" << std::endl;
        std::clog << "-- Name of port: " << get_name() << std::endl;
        std::clog << "-- Copy values to local buffer: " << copy_values << std::endl;
        std::clog << "-- Number of input elements: " << n_input << std::endl;
#endif
        if (is_fixed()) {
#if IMPBFF_VERBOSE
            std::clog << "WARNING: The port is fixed the action will be ignored." << std::endl;
#endif
            return;
        }
        update_value_type<T>(n_input);
#if IMPBFF_VERBOSE
        std::clog << "-- Copying values to local buffer." << std::endl;
#endif
        cc_.set(input, n_input);
        if (node_ != nullptr) {
#if IMPBFF_VERBOSE
            std::clog << "-- Updating attached node." << std::endl;
#endif
            update_attached_node();
#if IMPBFF_VERBOSE
            std::clog << "-- Updating dependent ports." << std::endl;
#endif
            set_value_of_dependents(input, n_input);
        }
    }

    template<typename T>
    void update_buffer() {
        auto v = get_array<T>("value");
        cc_.set(v.data(), v.size());
    }

    template<typename T>
    void get_own_value(T **output, int *n_output, bool update_local_buffer = false) {
        if ((cc_.size() == 0) || update_local_buffer) {
            update_buffer<T>();
        }
        T* origin = (T *) malloc(cc_.size());
        auto n_buffer_elements = cc_.size() / sizeof(T);
        memcpy(origin, cc_.get_bytes().data(), cc_.size());
        if (is_bounded() && bound_is_valid())
            Functions::bound_values<T>(
                    origin, n_buffer_elements, bounds_[0], bounds_[1]);
        *n_output = n_buffer_elements;
        *output = origin;
    }

    template<typename T>
    void get_value(T **output, int *n_output) {
        if (!is_linked()) {
#if IMPBFF_VERBOSE
            std::clog << "GET VALUE" << std::endl;
            std::clog << "-- Name of Port: " << get_name() << std::endl;
            std::clog << "-- Local value type: " << value_type << std::endl;
            std::clog << "-- Local buffer is filled: " << (n_buffer_elements_ > 0) << std::endl;
            std::clog << "-- Number of elements in local buffer: " << n_buffer_elements_ << std::endl;
            std::clog << "-- Port is not linked." << std::endl;
#endif
            get_own_value(output, n_output);
        } else {
#if IMPBFF_VERBOSE
            std::clog << "GET VALUE" << std::endl;
            std::clog << "-- Port is linked to " << get_link()->get_name() << std::endl;
#endif
            get_link()->get_value(output, n_output);
        }
#if IMPBFF_VERBOSE
        std::clog << "-- Number of elements: " << *n_output << std::endl;
#endif
    }

    virtual bson_t get_bson() final;

    // Methods
    //--------------------------------------------------------------------
    void update_attached_node();

    void set_bounded(bool bounded);

    bool bound_is_valid();

    bool is_bounded();

    void set_bounds(double *input, int n_input);

    void get_bounds(double **output, int *n_output);

    bool is_fixed();

    void set_fixed(bool fixed);

    bool is_output();

    void set_port_type(bool is_output);

    bool is_reactive();

    bool is_float();

    void set_reactive(bool reactive);

    bool write_to_db();

    bool read_from_db(const std::string &oid_string);

    void get_bytes(unsigned char **output, int *n_output, bool copy = false);

    void set_bytes(unsigned char *input, int n_input);
    
    void set_link(std::shared_ptr<CnPort> v);

    bool unlink();

    bool is_linked();

    std::vector<CnPort *> get_linked_ports();

    std::shared_ptr<CnPort> get_link();

};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_CNPORT_H
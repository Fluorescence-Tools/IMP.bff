#ifndef IMPBFF_PORT_H
#define IMPBFF_PORT_H

#include <IMP/bff/bff_config.h>
#include <IMP/object_macros.h>

#include <cstdint>
#include <memory>
#include <algorithm> /* std::max, std::min */
#include <cmath>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>

#include <IMP/bff/MongoObject.h>
#include <IMP/bff/Node.h>
#include <IMP/bff/Functions.h>

IMPBFF_BEGIN_NAMESPACE
class Node;

/// Specifies the type of the Port
/*!
 * 0: long vector
 * 1: double vector
 * 2: numpy binary
 * 3: long single number
 * 4: double single number
 */
typedef enum{
    BFF_PORT_VECTOR_INT,
    BFF_PORT_VECTOR_FLOAT,
    BFF_PORT_INT,
    BFF_PORT_FLOAT,
    BFF_PORT_BINARY
} PortType;

typedef enum{
    BFF_PORT_ADD,
    BFF_PORT_MUL
} PortOperation;


class IMPBFFEXPORT Port : public MongoObject{

private:

    void *buffer_ = nullptr;
    int n_buffer_elements_ = 0;
    int buffer_element_size_ = 1;
    std::vector<double> bounds_{};
    std::shared_ptr<Node> node_ = nullptr;
    bool own_buffer = false;

    /*!
     * @brief This attribute can point to another Port (default value nullptr).
     * If this is not a nullptr the value returned by the Port is the value
     * it the value the other Port.
     */
    std::shared_ptr<Port> link_ = nullptr;

    /*!
     * @brief This attribute stores the Ports that are dependent on the value
     * of this Port object. If this Port object is reactive a change of the
     * value of this Port object is propagated to the dependent Ports.
     */
    std::vector<Port *> linked_to_;

    bool remove_links_to_port();

    template<typename T>
    void set_value_of_dependents(T *input, int n_input) {
        for (auto &v : linked_to_) {
            v->set_value(input, n_input);
        }
    }

    int value_type = BFF_PORT_VECTOR_INT;

public:

    size_t size() {
        return n_buffer_elements_;
    }

    // Constructor & Destructor
    //--------------------------------------------------------------------
    ~Port() {
        remove_links_to_port();
        if (own_buffer) {
            free(buffer_);
        }
    }

    Port(
            bool fixed = false,
            bool is_output = false,
            bool is_reactive = false,
            bool is_bounded = false,
            double lb = 0,
            double ub = 0,
            int value_type = BFF_PORT_VECTOR_FLOAT,
            std::string name = ""
    ) : MongoObject(name) {
        append_string(&document, "type", "port");
        bson_append_bool(&document, "is_output", 9, false);
        bson_append_bool(&document, "is_fixed", 8, false);
        bson_append_bool(&document, "is_reactive", 11, false);
        bson_append_oid(&document, "link", 4, &oid_document);
        bson_append_bool(&document, "is_bounded", 10, false);
        buffer_ = malloc(1);

        set_fixed(fixed);
        set_port_type(is_output);
        set_reactive(is_reactive);
        set_bounded(is_bounded);
        if (is_bounded) {
            bounds_.push_back(lb);
            bounds_.push_back(ub);
        }
        Port::value_type = value_type;
    }

    void set_node(std::shared_ptr<Node> node_ptr);

    std::shared_ptr<Node> get_node();

    // Getter & Setter
    //--------------------------------------------------------------------
    int get_value_type() {
        if (!is_linked()) {
            return value_type;
        } else {
            return get_link()->value_type;
        }
    }

    void set_value_type(int v) {
        value_type = v;
        if (is_linked()) get_link()->set_value_type(v);
    }

    template<typename T>
    void set_value(T *input, int n_input, bool copy_values = true) {
#if IMPBFF_VERBOSE
        std::clog << "SET PORT VALUE" << std::endl;
        std::clog << "-- Name of port: " << get_name() << std::endl;
        std::clog << "-- Copy values to local buffer: " << copy_values << std::endl;
        std::clog << "-- Number of input elements: " << n_input << std::endl;
#endif
        if (is_fixed()) {
#if IMPBFF_VERBOSE
            std::clog << "WARNING: The port is fixed the action is ignored." << std::endl;
#endif
            return;
        }
        own_buffer = copy_values;
        if (std::is_same<T, double>::value) {
            value_type = 1;
        } else {
            value_type = 0;
        }
        // use scalar type for single numbers
        if(n_input == 1) value_type += 2;
        if(!copy_values){
#if IMPBFF_VERBOSE
            std::clog << "-- Avoiding copy - assign pointer of local buffer to input." << std::endl;
#endif
            // use pointer
            free(buffer_);
            buffer_ = input;
            buffer_element_size_ = sizeof(T);
            n_buffer_elements_ = n_input;
        } else{
#if IMPBFF_VERBOSE
            std::clog << "-- Copying values to local buffer." << std::endl;
#endif
            if ((int) (n_input * sizeof(T)) > (int) (n_buffer_elements_ * buffer_element_size_)){
#if IMPBFF_VERBOSE
                std::clog << "-- Size of input exceeds the local buffer: reallocating" << std::endl;
#endif
                buffer_ = std::realloc(buffer_, n_input * sizeof(T));
                buffer_element_size_ = sizeof(T);
            }
            if (buffer_ != nullptr) {
                memcpy(buffer_, input, n_input * sizeof(T));
            } else{
                std::cerr << "ERROR: Could not reallocate local buffer. " << std::endl;
            }
            n_buffer_elements_ = n_input;
        }
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
        n_buffer_elements_ = v.size();
        buffer_ = realloc(buffer_, buffer_element_size_ * n_buffer_elements_);
        memcpy(buffer_, v.data(), n_buffer_elements_ * sizeof(T));
    }

    template<typename T>
    void get_own_value(T **output, int *n_output, bool update_local_buffer = false) {
#if IMPBFF_VERBOSE
        std::clog << "GET OWN VALUE" << std::endl;
        std::clog << "-- Name of Port: " << get_name() << std::endl;
        std::clog << "-- Local value type: " << value_type << std::endl;
        std::clog << "-- Local buffer is filled: " << (n_buffer_elements_ > 0) << std::endl;
        std::clog << "-- Update local buffer: " << update_local_buffer << std::endl;
        std::clog << "-- Number of elements in local buffer: " << n_buffer_elements_ << std::endl;
#endif
        if ((n_buffer_elements_ == 0) || update_local_buffer) {
#if IMPBFF_VERBOSE
            std::clog << "-- Updating local buffer from bson" << std::endl;
#endif
            update_buffer<T>();
        }
#if IMPBFF_VERBOSE
        std::clog << "-- Checking if types are matching:" << std::endl;
#endif
        bool types_match = (
                ((std::is_same<T, double>::value) && ((value_type == 1) || (value_type == 3))) ||
                ((std::is_same<T, long>::value) && ((value_type == 0) || (value_type == 2)))
        );
        T* origin;
        if(!types_match){
#if IMPBFF_VERBOSE
            std::clog << "-- Requested type does not match buffer type - recasting types." << std::endl;
#endif
            origin = (T *) malloc(n_buffer_elements_ * sizeof(T));
            if((value_type == 1) || (value_type == 3)){
                auto b = reinterpret_cast<double *>(buffer_);
                for(int i = 0; i < n_buffer_elements_; i++){
                    origin[i] = (T) b[i];
                }
            }
            if((value_type == 0) || (value_type == 2)){
                auto b = reinterpret_cast<long *>(buffer_);
                for(int i = 0; i < n_buffer_elements_; i++){
                    origin[i] = (T) b[i];
                }
            }
        } else{
#if IMPBFF_VERBOSE
            std::clog << "-- Requested type matches local type. " << std::endl;
#endif
            origin = reinterpret_cast<T *>(buffer_);
        }
#if IMPBFF_VERBOSE
        std::clog << "-- Checking if value is bounded:" << std::endl;
#endif
        if (is_bounded() && bound_is_valid()) {
#if IMPBFF_VERBOSE
            std::clog << "-- Values are bounded " << std::endl;
            std::clog << "-- Bounds: [" << bounds_[0] << ", " << bounds_[1] << "]" << std::endl;
#endif
            auto bounded_array = (T *) malloc(n_buffer_elements_ * sizeof(T));
            memcpy(bounded_array, origin, (size_t) n_buffer_elements_ * sizeof(T));
            Functions::bound_values<T>(bounded_array, n_buffer_elements_, bounds_[0], bounds_[1]);
            *n_output = n_buffer_elements_;
            *output = reinterpret_cast<T *>(bounded_array);
        } else {
#if IMPBFF_VERBOSE
            std::clog << "-- Values are not bounded " << std::endl;
#endif
            *n_output = n_buffer_elements_;
            *output = origin;
        }
    }

    template<typename T>
    void get_value(T **output, int *n_output) {
        if (!is_linked()) {
#if IMPBFF_VERBOSE
            std::clog << "GET VALUE" << std::endl;
            std::clog << "-- Name of Port: " << get_name() << std::endl;
            std::clog << "-- Local value type: " << value_type << std::endl;
            std::clog << "-- Local buffer filled: " << (n_buffer_elements_ > 0) << std::endl;
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

    void set_bounds(std::vector<double> bounds);

    std::vector<double> get_bounds();

    bool is_fixed();

    void set_fixed(bool fixed);

    bool is_output();

    void set_port_type(bool is_output);

    bool is_reactive();

    bool is_float() {
        return ((get_value_type() == BFF_PORT_VECTOR_FLOAT) || (get_value_type() == BFF_PORT_FLOAT));
    }

    void set_reactive(bool reactive);

    bool write_to_db();

    bool read_from_db(const std::string &oid_string);

    std::vector<unsigned char> get_bytes();

    void set_bytes(unsigned char *input, int n_input);

    void set_buffer_ptr(size_t ptr, int n_elements, int element_size) {
        buffer_ = (void *) ptr;
        buffer_element_size_ = element_size;
        n_buffer_elements_ = n_elements;
    }

    size_t get_buffer_ptr() {
        return (size_t) (&buffer_);
    }

    void set_link(std::shared_ptr<Port> v) {
#if IMPBFF_VERBOSE
        std::clog << "SET_LINK" << std::endl;
        std::clog << "-- Link to Port: " << v->get_name() << std::endl;
        std::clog << "-- Link value type: " << v->value_type << std::endl;
#endif
        if (v != nullptr) {
            if((v.get() == this)){
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

    bool unlink() {
        set_oid("link", get_bson_oid());
        bool re = true;
        re &= remove_links_to_port();
        link_ = nullptr;
        return re;
    }

    bool is_linked() {
        return (link_ != nullptr);
    }

    std::vector<Port *> get_linked_ports() {
        return linked_to_;
    }

    std::shared_ptr<Port> get_link() {
        return link_;
    }


};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PORT_H
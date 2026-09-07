/*
    A minimal map-like container that preserves insertion order, for use
    as nlohmann::basic_json's ObjectType.

    Backport of nlohmann::ordered_map (JSON for Modern C++, Niels Lohmann,
    MIT License, https://github.com/nlohmann/json) into the vendored 3.6.1
    single-header in include/internal/json.h, which predates the type. The
    class body is upstream's, trimmed to the plain key_type overloads the
    3.6.1 machinery calls; the SFINAE-keyed heterogeneous-lookup overloads
    of newer upstreams are left out because 3.6.1 has no
    is_usable_as_key_type to guard them with.

    Session.cpp needs it because chinet's session format is order-carrying:
    a node's "ports" map restores the node's port insertion order, and that
    order is observable through the operator callbacks (they compute over
    the first two inputs in insertion order). nlohmann's default object
    type is std::map, which sorts keys and so scrambles the document order
    on load -- json.loads on the Python side of the A/B tests keeps it, and
    the ported runtime must too.
*/

#ifndef IMPBFF_INTERNAL_ORDERED_MAP_H
#define IMPBFF_INTERNAL_ORDERED_MAP_H

#include <algorithm>  // std::equal
#include <cstddef>  // size_t
#include <functional>  // less, equal_to
#include <initializer_list>  // initializer_list
#include <iterator>  // input_iterator_tag, iterator_traits
#include <memory>  // allocator
#include <stdexcept>  // out_of_range
#include <type_traits>  // enable_if, is_convertible
#include <utility>  // pair, move, forward
#include <vector>  // vector

namespace nlohmann {

template <class Key, class T, class IgnoredLess = std::less<Key>,
          class Allocator = std::allocator<std::pair<const Key, T> > >
struct ordered_map
    : std::vector<std::pair<const Key, T>, Allocator> {
  using key_type = Key;
  using mapped_type = T;
  using Container = std::vector<std::pair<const Key, T>, Allocator>;
  using iterator = typename Container::iterator;
  using const_iterator = typename Container::const_iterator;
  using size_type = typename Container::size_type;
  using value_type = typename Container::value_type;
  using key_compare = std::equal_to<Key>;

  // Explicit constructors instead of `using Container::Container` --
  // older compilers (GCC <= 5.5, Xcode <= 9.4) choke on the inherited ones.
  ordered_map() noexcept(noexcept(Container())) : Container() {}
  explicit ordered_map(const Allocator& alloc) noexcept(
      noexcept(Container(alloc)))
      : Container(alloc) {}
  template <class It>
  ordered_map(It first, It last, const Allocator& alloc = Allocator())
      : Container(first, last, alloc) {}
  ordered_map(std::initializer_list<value_type> init,
              const Allocator& alloc = Allocator())
      : Container(init, alloc) {}

  std::pair<iterator, bool> emplace(const key_type& key, T&& t) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return std::make_pair(it, false);
      }
    }
    Container::emplace_back(key, std::forward<T>(t));
    return std::make_pair(std::prev(this->end()), true);
  }

  T& operator[](const key_type& key) {
    return emplace(key, T()).first->second;
  }

  const T& operator[](const key_type& key) const { return at(key); }

  T& at(const key_type& key) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return it->second;
      }
    }
    throw std::out_of_range("key not found");
  }

  const T& at(const key_type& key) const {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return it->second;
      }
    }
    throw std::out_of_range("key not found");
  }

  size_type erase(const key_type& key) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        // const Keys cannot move, so re-construct them in place.
        for (auto next = it; ++next != this->end(); ++it) {
          it->~value_type();  // destroy but keep the allocation
          new (&*it) value_type(std::move(*next));
        }
        Container::pop_back();
        return 1;
      }
    }
    return 0;
  }

  iterator erase(iterator pos) { return erase(pos, std::next(pos)); }

  iterator erase(iterator first, iterator last) {
    if (first == last) {
      return first;
    }
    const typename Container::difference_type elements_affected =
        std::distance(first, last);
    const typename Container::difference_type offset =
        std::distance(Container::begin(), first);
    for (auto it = first;
         std::next(it, elements_affected) != Container::end(); ++it) {
      it->~value_type();  // destroy but keep the allocation
      new (&*it) value_type(
          std::move(*std::next(it, elements_affected)));
    }
    Container::resize(this->size() - static_cast<size_type>(elements_affected));
    return Container::begin() + offset;
  }

  size_type count(const key_type& key) const {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return 1;
      }
    }
    return 0;
  }

  iterator find(const key_type& key) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return it;
      }
    }
    return Container::end();
  }

  const_iterator find(const key_type& key) const {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, key)) {
        return it;
      }
    }
    return Container::end();
  }

  std::pair<iterator, bool> insert(value_type&& value) {
    return emplace(value.first, std::move(value.second));
  }

  std::pair<iterator, bool> insert(const value_type& value) {
    for (auto it = this->begin(); it != this->end(); ++it) {
      if (m_compare(it->first, value.first)) {
        return std::make_pair(it, false);
      }
    }
    Container::push_back(value);
    return std::make_pair(--this->end(), true);
  }

  template <typename InputIt>
  using require_input_iter =
      typename std::enable_if<std::is_convertible<
          typename std::iterator_traits<InputIt>::iterator_category,
          std::input_iterator_tag>::value>::type;

  template <typename InputIt, typename = require_input_iter<InputIt> >
  void insert(InputIt first, InputIt last) {
    for (auto it = first; it != last; ++it) {
      insert(*it);
    }
  }

 private:
  key_compare m_compare;
};

}  // namespace nlohmann

#endif  // IMPBFF_INTERNAL_ORDERED_MAP_H

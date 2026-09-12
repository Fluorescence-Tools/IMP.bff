/**
 *  \file IMP/bff/GraphSession.h
 *  \brief Save and load node graphs in chinet's on-disk session format.
 *
 *  A GraphSession is what a chisurf .csp project is: a named graph of nodes
 *  (each with named ports) plus the free-floating ports chisurf's fit
 *  parameters are. save() writes it as chinet wrote it -- JSONL, one JSON
 *  document per line, the session line first, then every object's
 *  document, cross-referenced by _id -- and load() reads that format back,
 *  resolving session.nodes, node.ports and port.link by _id. The format is
 *  the spec: a chinet-generated fixture is committed at
 *  test/session/chinet_fixture.jsonl, and load() reconstructs it field for
 *  field. Field order within a line matches chinet's writes (the order of
 *  the lines across the file does not matter to a reader and only needs to
 *  stay deterministic, which it is: nodes in the order added, each node's
 *  ports after it, free ports last).
 *
 *  Ported from chinet's chinet/session.py (phase 2 of removing chinet from
 *  chisurf, on top of the phase-1 GraphPort/GraphNode runtime). Deliberately
 *  standalone, like the runtime: no IMP particles, restraints or
 *  decorators take part. What is not ported from chinet:
 *
 *  - the DB registry (db.py) and the MMFDB backend. chinet registered
 *    every constructed object with a global database and save() dumped
 *    whatever had accumulated there; here the session is the registry --
 *    a node or port is persisted only by add_node()/add_port(), never by
 *    construction alone, so no hidden global state crosses the wrapper.
 *    (Loading does attach every document it finds: unclaimed ports become
 *    the session's free ports, which is what chinet's DB held for them.)
 *  - chinet's write_to_db/read_from_db/connect_to_db object plumbing,
 *    append_object, the schema.py chinet.session.v1 conversion, and the
 *    create_port/create_node template helpers.
 *  - the legacy monolithic {"session": ..., "objects": [...]} format IS
 *    read (old projects on disk), only not written.
 *
 *  Semantics carried over exactly, because round-tripping a real file
 *  depends on them:
 *
 *  - loading restores a port's value, then its value_type: the type code
 *    is the document's, not inferred from the data (a saved int port with
 *    value 2 loads back as an int port, code 0, value 2).
 *  - a port's saved value follows its link, as chinet's .value does; the
 *    link field is the target port's _id.
 *  - int-typed ports (codes 0 and 2) save their values as JSON integers,
 *    float-typed ports as JSON floats -- chinet's documents carry numpy's
 *    int64/float64 distinction that way.
 *  - a node's saved "valid" is the raw flag, not is_valid(); restoring
 *    links in document order can invalidate a node again, exactly as
 *    chinet's load does.
 *  - _id and precursor round-trip verbatim: save-load-save reproduces the
 *    same document per _id.
 *
 *  Two small divergences, both noted where they happen: enabling bounds
 *  on load clips the loaded value into them (chinet's set_document did
 *  not clip; a chinet-written file is always within bounds, so this only
 *  bites hand-edited files), and a prior that is not JSON object text is
 *  refused at save rather than silently written as null.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_GRAPHSESSION_H
#define IMPBFF_GRAPHSESSION_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A saved-and-loadable graph of nodes and free ports, chinet's GraphSession.
/*!
    Sessions hold shared_ptr-owned nodes and ports; nothing is registered
    anywhere else, so a session that goes away takes its graph with it
    unless Python (or a link) still holds a reference.
*/
class IMPBFFEXPORT GraphSession : public GraphObject {
 public:
  //! Nodes by key (see get_nodes(); insertion order is kept internally).
  typedef std::map<std::string, std::shared_ptr<GraphNode> > GraphNodeMap;
  //! The free-floating ports, in the order added.
  typedef std::vector<std::shared_ptr<GraphPort> > GraphPortList;

  //! An empty session named "session", with a fresh uid.
  GraphSession();
  virtual ~GraphSession();

  //! Hold a node under a key (the session line's nodes map). Replaces a
  //! node already under that key.
  void add_node(const std::string& key, std::shared_ptr<GraphNode> node);
  //! All nodes by key.
  GraphNodeMap get_nodes() const;
  //! The node under a key, or a null pointer.
  std::shared_ptr<GraphNode> get_node(const std::string& key) const;

  //! Hold a free-floating port (chisurf's unattached fit parameters).
  /*!
      Ports of a held node are persisted with the node and need no
      add_port; this is for ports that belong to no node. Adding a port
      twice is a no-op.
  */
  void add_port(std::shared_ptr<GraphPort> port);
  //! The free-floating ports, in the order added.
  GraphPortList get_ports() const;
  //! The first free port with this name, or a null pointer.
  std::shared_ptr<GraphPort> get_port(const std::string& name) const;

  //! Let go of every node and free port (the counterpart of add_node /
  //! add_port; the objects themselves stay alive while anything holds
  //! them). The session's own uid is untouched.
  void clear();

  //! How many nodes the session holds.
  unsigned int get_number_of_nodes() const;
  //! How many distinct ports the session would write: every port of
  //! every node plus the free ports.
  unsigned int get_number_of_ports() const;

  //! Write the session and every object it holds to a JSONL file.
  /*!
      The session line first, then each node followed by its ports (in
      the order they were added), then the free ports; one document per
      _id. Throws std::runtime_error if the file cannot be written, and
      std::invalid_argument if a port's prior is not JSON object text.
  */
  void save(const std::string& path) const;

  //! Read a session back from a JSONL (or legacy monolithic) file.
  /*!
      Returns a null pointer for an empty file or a JSONL file without a
      session line, chinet's None; throws std::runtime_error for a file
      that is neither format or holds a malformed line. Documents are
      resolved by _id: session.nodes to nodes, node.ports to ports,
      port.link to its target port. Ports no node claims become the
      session's free ports, and uid/precursor are preserved exactly.
  */
  static std::shared_ptr<GraphSession> load(const std::string& path);

 private:
  GraphNodeMap nodes_;
  //! Keys of nodes_ in add order (std::map sorts; the session line's
  //! nodes map is written in the order the nodes were added).
  std::vector<std::string> node_order_;
  GraphPortList ports_;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_GRAPHSESSION_H

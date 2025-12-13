//
// Created by Mert Biyikli on 20.05.25.
//

#ifndef OBLIVIOUSROUTING_TREE_H
#define OBLIVIOUSROUTING_TREE_H

#include "Node.h"
#include <fstream>
#include <sstream>
#include <queue>
#include "../../raecke_oracle_iteration.h"
/*
class FRT_Tree {
    std::shared_ptr<FRT_Node> root;
    std::map<int, std::vector<std::shared_ptr<FRT_Node>>> idLayer2nodes;
    std::vector<std::shared_ptr<FRT_Node>> all_nodes;
public:
    void createRoot(int layers) {
        root = std::make_shared<FRT_Node>();

        // If no nodes in the layer, initialize a new vector for the layer
        if (idLayer2nodes[layers].empty()) {
            idLayer2nodes[layers] = {root};
        } else {
            idLayer2nodes[layers].push_back(root);  // Replace with root node if already exists
        }

    }

    std::shared_ptr<FRT_Node> GetRoot() {
        return root;
    }

    void SetRoot(std::shared_ptr<FRT_Node> r) {
        root = std::move(r);
    }

    void AddNodeToLayer(int layer, std::shared_ptr<FRT_Node> node) {
        if(node) {
            idLayer2nodes[layer].push_back(node);
        }
    }

    std::vector<std::shared_ptr<FRT_Node>>& GetLayerNodes(int layer) {
        return idLayer2nodes[layer];  // Return empty vector if layer not found
    }

    const auto& GetLayer2Nodes() const {
        return idLayer2nodes;
    }

    void print() {
        return;
        // print the tree in a BFS fashion starting from the root
        if (!root) {
            std::cout << "Tree is empty." << std::endl;
            return;
        }

        std::queue<std::shared_ptr<FRT_Node>> q;
        q.push(root);
        while (!q.empty()) {
            auto node = q.front();
            q.pop();
            // also print the layer of the nodes
            if (!node) continue;  // Skip null nodes
            // Print node information

            for(auto& [layer, nodes] : idLayer2nodes) {
                if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) {
                    std::cout << "Node at layer " << layer << ", ";
                    break;
                }
            }
            std::cout << "Node with Center: " << node->getCenter() << ", Vertices: ";
            for (const auto& v : node->GetVertices()) {
                std::cout << v << " ";
            }
            std::cout << std::endl;

            // TODO: FIX THIS
            for (const auto& child : node->getChildren()) {
                //q.push(dynamic_cast<FRT_Node>(child));
            }
        }
    }

    std::shared_ptr<FRT_Node> getNode(int id) const {
        if (id < 0 || id >= (int)all_nodes.size())
            throw std::runtime_error("Invalid node id");
        return all_nodes[id];
    }

    // --- Additions ---
    int addClusterNode(const std::vector<int>& vertices, int level = 0) {
        auto node = std::make_shared<FRT_Node>();

        for (auto& v : vertices) {
            node->addVertex(v);
        }
        if (!vertices.empty()) node->setCenter(vertices.front());

        idLayer2nodes[level].push_back(node);
        all_nodes.push_back(node);

        // --- Automatically set the first node as root ---
        if (!root) {
            root = node;
        }

        return (int)all_nodes.size() - 1;
    }

    void addTreeEdge(int child_id, int parent_id) {
        if (child_id < 0 || parent_id < 0 ||
            child_id >= (int)all_nodes.size() || parent_id >= (int)all_nodes.size())
            throw std::runtime_error("Invalid node id in addTreeEdge.");

        auto parent = all_nodes[parent_id];
        auto child = all_nodes[child_id];

        child->setParent(parent);
        parent->addChild(child);
    }
};

class Tree {
public:
    std::shared_ptr<Node> root;
    int id;  // Unused, but can be kept for consistency
    std::map<int, std::vector<std::shared_ptr<Node>>> idLayer2nodes;

    Tree() : root(nullptr), id(0) {}

    void clear() {
        root.reset();  // Clear the root node
        idLayer2nodes.clear();  // Clear the layer to nodes mapping
        id = 0;  // Reset the ID
    }

    std::shared_ptr<Node> GetRoot() {
        return root;
    }

    void SetRoot(std::shared_ptr<Node> root) {
        this->root = std::move(root);
    }

    int GetLayerCount() const {
        return idLayer2nodes.size();
    }

    int GetId() const {
        return id;
    }

    void SetId(int id) {
        this->id = id;
    }

    void CreateRoot(int layers) {
        this->root = std::make_shared<Node>();  // Create root as shared_ptr
        Node::nextId = 0;  // Reset the static ID counter for Node
        this->root->SetId(Node::nextId);  // Set the ID of the root node to 0
        Node::nextId++;

        // If no nodes in the layer, initialize a new vector for the layer
        if (idLayer2nodes[layers].empty()) {
            idLayer2nodes[layers] = {root};
        } else {
            idLayer2nodes[layers].push_back(root);  // Replace with root node if already exists
        }

    }

    bool IsLayerWithAllSingletons(int i) const {
        bool isAllSingletons = true;
        for (const auto& node : idLayer2nodes.at(i)) {
            if (node->GetVertices().size() != 1) {
                isAllSingletons = false;
                break;
            }
        }
        return isAllSingletons;
    }


    // Get all nodes in a layer
    std::vector<std::shared_ptr<Node>> GetLayerNodes(int layer) {
        auto it = idLayer2nodes.find(layer);
        return it != idLayer2nodes.end() ? it->second : std::vector<std::shared_ptr<Node>>();  // Return empty vector if layer not found
    }

    std::vector<std::shared_ptr<Node>> GetLayerNodes(int layer) const{
        auto it = idLayer2nodes.find(layer);
        return it != idLayer2nodes.end() ? it->second : std::vector<std::shared_ptr<Node>>();  // Return empty vector if layer not found
    }



    void AddNodeToLayer(int i, std::shared_ptr<Node> childNode) {
        if(idLayer2nodes[i].empty()) {
            this->idLayer2nodes[i] = {childNode};
        }else{
            idLayer2nodes[i].push_back(childNode);
        }

    }

    std::string toString() const {
        return root ? root->toString(0) : "No root";
    }

    std::shared_ptr<Node> GetLeafVertex(int id) const {
        for (const auto& node : idLayer2nodes.at(0)) {
            if (node->ContainsVertex(id)) {
                return node;
            }
        }
        return nullptr;
    }

    bool operator<(const Tree& other) const {
        return this->id > other.id;
    }

    void ExportToDot(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return;

        file << "digraph Tree {\n";
        file << "  node [shape=ellipse, fontsize=12];\n";

        for (const auto& [layer, nodes] : idLayer2nodes) {
            for (const auto& node : nodes) {
                std::stringstream nodeLabel;
                nodeLabel << "Layer " << layer << "\\nID: " << node->GetId();
                nodeLabel << "\\nCenter: " << node->GetCenter();
                nodeLabel << "\\nVertices: ";
                for (int v : node->GetVertices()) nodeLabel << v << ",";
                file << "  \"" << node->GetId() << "\" [label=\"" << nodeLabel.str() << "\"];\n";

                if (auto parent = node->GetParent()) {
                    file << "  \"" << parent->GetId() << "\" -> \"" << node->GetId() << "\";\n";
                }
            }
        }

        file << "}\n";
        file.close();
    }

};*/



#endif //OBLIVIOUSROUTING_TREE_H

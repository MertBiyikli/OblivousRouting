//
// Created by Mert Biyikli on 20.05.25.
//

#ifndef OBLIVIOUSROUTING_NODE_H
#define OBLIVIOUSROUTING_NODE_H

#include <memory>
#include <vector>
#include <set>
#include <algorithm>

#include "../../raecke_oracle_iteration.h"

class EfficientFRTTreeNode : public ITreeNode {
public:
    int center = -1;
    int id = -1;
    std::vector<int> members; // original graph nodes
    std::weak_ptr<EfficientFRTTreeNode> parent;
    std::vector<std::shared_ptr<EfficientFRTTreeNode>> children;

    virtual std::shared_ptr<ITreeNode> getParent() const override {
        return parent.lock();
    }

    std::vector<std::shared_ptr<ITreeNode>> getChildren() override {
        std::vector<std::shared_ptr<ITreeNode>> child_ptrs;
        for (std::shared_ptr<EfficientFRTTreeNode>& child : children) {
            child_ptrs.push_back(child);
        }
        return child_ptrs;
    }

    const std::vector<int> &getMembers() const override {
        return members;
    }

};
/*
class FRT_Node :public ITreeNode{
public:
    int id;
    std::weak_ptr<FRT_Node> parent;
    std::vector<std::shared_ptr<FRT_Node>> children_;
    std::set<int> vertices;
    int center;

    FRT_Node() : id(-1), center(-1) {}

    ITreeNode* getParent() const override {
        return parent.lock().get();
    }
    const std::vector<int>& getMembers() const override {
        return std::vector<int>(vertices.begin(), vertices.end());
    }

     std::vector<ITreeNode*> getChildren() override {
        std::vector<ITreeNode*> child_ptrs;
        for (std::shared_ptr<FRT_Node>& child : children_) {
            child_ptrs.push_back(child.get());
        }
        return child_ptrs;
    }

    void addVertex(int vertex) {
        vertices.insert(vertex);
    }

    void addVertices(const std::vector<int>& v) {
        vertices.insert(v.begin(), v.end());
    }

    bool containsVertex(int vertex) const {
        return vertices.find(vertex) != vertices.end();
    }

    std::set<int> GetVertices() {
        return vertices;
    }

    void setCenter(int vertex) {
        center = vertex;
    }

    int getCenter() const {
        return center;
    }

    void setParent(std::shared_ptr<FRT_Node> p) {
        parent = p;
    }

    void addChild(std::shared_ptr<FRT_Node> child) {
        children_.push_back(child);
    }

    const std::vector<std::shared_ptr<FRT_Node>>& getChildren() const {
        return children_;
    }

    std::string toString(int levels = 0) const {
        std::string result = std::string(levels, '.');
        result += "Center: " + std::to_string(center) + " Vertices: ";
        for (const auto& v : vertices) {
            result += std::to_string(v) + " ";
        }
        result += "\n";

        for (const auto& child : children_) {
            result += child->toString(levels + 1);
        }

        return result;
    }


};

class Node {
public:
        static int nextId;
        int id;
        std::weak_ptr<Node> parent;
        std::vector<std::shared_ptr<Node>> children;
        std::vector<int> vertices;
        int center;

        Node() :
        id(nextId++), center(-1){
        }

        Node(const Node& other) :
        id(other.id), center(other.center), parent(other.parent),
        children(other.children), vertices(other.vertices) {
            // Copy constructor
        }

        // Getter and Setter for center
        int GetCenter() const{
            return center;
        }

        void SetCenter(int c) {
            center = c;
        }

        // Getter and Setter for id
        int GetId() const {
            return id;
        }

        void SetId(int id) {
            this->id = id;
        }

        // Getter and Setter for parent
        std::shared_ptr<Node> GetParent() const {
            return parent.lock();
        }

        void SetParent(std::shared_ptr<Node> p) {
            parent = p;
        }

        // Getter and Setter for children
        const std::vector<std::shared_ptr<Node>>& GetChildren() const {
            return children;
        }

        void SetChildren(std::vector<std::unique_ptr<Node>>& c) {
            for(auto& node : c) {
                children.push_back(std::move(node));
            }
        }

        // Getter and Setter for vertices
        const std::vector<int>& GetVertices() const {
            return vertices;
        }

        void SetVertices(const std::vector<int>& v) {
            vertices = v;
        }

        void AddVertex(const int& vertex) {
            vertices.push_back(vertex);
        }

        void AddVertices(const std::vector<int>& v) {
            vertices.insert(vertices.end(), v.begin(), v.end());
        }

        bool ContainsVertex(int id) const {
            auto it = std::find_if(vertices.begin(), vertices.end(), [id](const int& v) {
                return v == id;
            });

            return it != vertices.end();
        }


        void AddChild(std::shared_ptr<Node> child) {
            children.push_back(child);
        }

        std::string toString(int levels = 0) const {
            std::string result = std::string(levels, '.');
            result += "Center: " + std::to_string(center ? center : -1) + " Vertices: ";
            for (const auto& v : vertices) {
                result += std::to_string(v) + " ";
            }
            result += "\n";

            for (const auto& child : children) {
                result += child->toString(levels + 1);
            }

            return result;
        }

        std::string toString() const {
            std::string result = "Center: " + std::to_string(center ? center : -1) + " Vertices: ";
            for (const auto& v : vertices) {
                result += std::to_string(v) + " ";
            }
            return result;
        }
};*/


#endif //OBLIVIOUSROUTING_NODE_H

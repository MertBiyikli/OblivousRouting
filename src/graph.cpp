//
// Created by Mert Biyikli on 11.05.25.
//


#include "graph.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <limits>


void Graph::readGraph(const std::string &filename) {
    std::ifstream infile(filename);
    std::cout << "Reading graph from file: " << filename << std::endl;
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int n = 0, m = 0;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'p') {
            std::string tmp;
            iss >> tmp >> tmp >> n >> m;
            break;
        }
    }

    if (n == 0) {
        throw std::runtime_error("Invalid or missing problem line in DIMACS file");
    }


    this->m_adj.resize(n);
    this->m_adj_capacities.resize(n);
    this->m_adj_distances.resize(n);

    do {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'a') {
            char dummy;
            int u, v;
            double cap = 1.0;
            iss >> dummy >> u >> v;
            if (!(iss >> cap)) cap = 1.0;
            this->addEdge(u - 1, v - 1, cap);  // convert to 0-based
        }
    } while (std::getline(infile, line));
}

void Graph::readLFGFile(const std::string &filename, bool withDistances) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open LGF file: " + filename);
    }

    std::string line;
    bool inNodesSection = false;
    bool inArcsSection = false;

    int maxNodeIdSeen = -1;
    int costColIdx = -1;
    int capColIdx  = -1;

    // Phase 1: Read lines until we finish @nodes and find @arcs header.
    std::vector<std::string> allLines;
    while (std::getline(infile, line)) {
        // Trim leading/trailing whitespace:
        // (Here we just check if line is empty or starts with a comment character)
        if (line.empty()) {
            continue;
        }
        allLines.push_back(line);
    }
    infile.close();

    // 1) First pass: scan for "@nodes" → read node lines → track maxNodeIdSeen.
    for (size_t i = 0; i < allLines.size(); ++i) {
        const std::string &raw = allLines[i];
        std::string lower = raw;
        // normalize to lowercase to match "@nodes" or "@arcs"
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (!inNodesSection) {
            if (lower.rfind("@nodes", 0) == 0) {
                inNodesSection = true;
                continue;
            }
        } else {
            // We are inside the @nodes block. Stop if we see "@arcs".
            if (lower.rfind("@edges", 0) == 0) {
                inNodesSection = false;
                inArcsSection = true;
                // The next line after this is expected to be the arcs‐header.
                continue;
            }
            // Otherwise, we expect this line to have: "<label><TAB><node_id>"
            // We split on tabs:
            std::istringstream iss(raw);
            std::vector<std::string> tokens;
            std::string tok;
            while (iss >> tok) {
                tokens.push_back(tok);
            }
            if (tokens.size() >= 2) {
                // node_id is the second column
                try {
                    int nid = std::stoi(tokens[1]);
                    maxNodeIdSeen = std::max(maxNodeIdSeen, nid);
                } catch (...) {
                    // malformed line: skip
                }

                try {
                    int nid = std::stoi(tokens[0]);
                    maxNodeIdSeen = std::max(maxNodeIdSeen, nid);
                } catch (...) {
                    // malformed line: skip
                }
            }
        }
    }

    if (maxNodeIdSeen < 0) {
        throw std::runtime_error("LGF file had no valid @nodes section or no node IDs found");
    }

    // 2) Initialize our graph with (maxNodeIdSeen+1) nodes
    int n = maxNodeIdSeen + 1;
    int m = 0;
    m_adj.clear();
    m_adj.resize(n);
    m_adj_distances.resize(n);
    m_adj_capacities.resize(n);
    m_vertices.resize(n);
    std::iota(m_vertices.begin(), m_vertices.end(), 0); // fill vertices with 0, 1, ..., maxNodeIdSeen

    // 3) Second pass: scan for "@arcs", find header, then parse each arc line:
    inArcsSection = false;
    bool readHeader = false;
    for (size_t i = 0; i < allLines.size(); ++i) {
        const std::string &raw = allLines[i];
        std::string lower = raw;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (!inArcsSection) {
            // Look for “@arcs”
            if (lower.rfind("@edges", 0) == 0) {
                inArcsSection = true;
            }
            continue;
        }

        // We are in @arcs section. The very first non‐empty line (after "@arcs") is the header.
        if (inArcsSection && !readHeader) {
            // This is the header line, e.g. "    label    cost    capacity"
            std::istringstream iss(raw);
            std::vector<std::string> headerTokens;
            std::string tok;
            while (iss >> tok) {
                if(!tok.empty()) {
                    headerTokens.push_back(tok);
                }
            }
            // Find which index is “cost” and which is “capacity”
            for (int col = 0; col < (int)headerTokens.size(); ++col) {
                std::string h = headerTokens[col];
                std::transform(h.begin(), h.end(), h.begin(),
                               [](unsigned char c) { return std::tolower(c); });
                if (h == "cost") {
                    costColIdx = col;
                }
                if (h == "capacity") {
                    capColIdx = col;
                }
            }
            costColIdx += 2; // TODO: this is hardcoded for now, should be dynamic
            readHeader = true;
            continue;
        }

        // Now each subsequent line (after header) is “u<TAB>v<TAB>…”
        if (raw.empty() || raw[0] == '#') {
            // skip blank or comment lines
            continue;
        }

        // Split on tabs
        std::istringstream iss(raw);
        std::vector<std::string> parts;
        std::string field;

        while (iss >> field) {
            if(!field.empty()) {
                parts.push_back(field);
            }
        }

        if (parts.size() < 2) {
            // not enough columns ⇒ skip
            continue;
        }

        // Parse u, v (1‐based IDs in LGF)
        int u, v;
        try {
            u = std::stoi(parts[0]);  // convert to 0‐based
            v = std::stoi(parts[1]);
        } catch (...) {
            continue; // malformed
        }
        if (u < 0 || v < 0 || u >= n || v >= n) {
            // invalid node Id ⇒ skip
            continue;
        }
        // We only add each undirected edge once (u < v)
        if (u > v) {
            continue;
        }

        // Determine capacity:
        double capacityValue = 1.0;
        bool capacityParsed = false;
        if (costColIdx >= 0 && costColIdx < (int)parts.size()) {
            // attempt to parse capacity
            try {
                capacityValue = std::stod(parts[costColIdx]);
                capacityParsed = true;
            } catch (...) {
                capacityParsed = false;
            }
        }
        if (!capacityParsed && withDistances && costColIdx >= 0 && costColIdx < (int)parts.size()) {
            // parse cost, then invert it
            try {
                double costVal = std::stod(parts[costColIdx]);
                if (costVal != 0.0) {
                    capacityValue = costVal;
                } else {
                    capacityValue = 1.0;
                }
                capacityParsed = true;
            } catch (...) {
                capacityParsed = false;
            }
        }
        // if still not parsed, we leave capacityValue = 1.0 by default

        // Finally, add the undirected edge (Graph::addEdge adds both directions)
        this->addEdge(u, v, capacityValue);
    }

    // 4) (Optional) If you have a PruneGraph class, prune degree‐one vertices
    //     new PruneGraph().pruneDegreeOne(*this);
    //
    // If you do not have PruneGraph or you don’t want to prune, simply comment out the line above.
}


/*
int Graph::numNodes() const {
    return m_iNumNodes;
}

int Graph::numEdges() const {
    return m_iNumEdges;
}

const std::vector<int>& Graph::GetVertices() const{
    return m_vertices;
}

std::vector<int>& Graph::GetVertices() {
    return m_vertices;
}


void Graph::addEdge(int source, int target, double capacity) {
    if (source < 0 || source >= m_iNumNodes || target < 0 || target >= m_iNumNodes) {
        throw std::out_of_range("Invalid node index");
    }

    auto edge = std::make_shared<Edge>();
    edge->source = source;
    edge->target = target;
    edge->capacity = capacity;

    auto rev_edge = std::make_shared<Edge>();
    rev_edge->source = target;
    rev_edge->target = source;
    rev_edge->capacity = capacity;
    m_adj[source].push_back(edge);
    m_adj[target].push_back(rev_edge);
    m_iNumEdges++;
}

*/



// … existing methods (numNodes, numEdges, addEdge, print, readGraph, GetDiGraph, etc.) …


// TODO: make the lgf reader more robust against variying column orderings
//       as well as more flexible with respect to the header lines.(it should accept @edges and @arcs interchangeably)

/**
 * Reads a graph in “.lgf” format.
 *
 * We expect the file to look roughly like:
 *
 * @nodes
 *    label    node_id
 *    foo      1
 *    bar      2
 *    baz      3
 * @arcs
 *        label    cost    capacity
 *    1    2       A1      5.0     10.0
 *    1    3       A2      2.5     8.0
 *    2    3       A3      7.0     12.0
 *
 * In other words:
 *  - “@nodes” section: each line has at least “<label><TAB><node_id>”. We ignore <label>, but track the largest node_id.
 *  - “@arcs” section: the first non‐empty, non‐comment line after “@arcs” is a header. We detect which column index is “cost” and which is “capacity”.
 *    Then each subsequent line has at least “<u><TAB><v>…”. We parse u,v (1‐based in typical LGF) → convert to 0‐based.
 *    If (u > v), we skip that line (because we only want each undirected edge once).
 *    If capacity column is present, we parse capacity. Otherwise if withDistances==true and cost column present, we do capacity = 1.0/cost.
 *    Otherwise default capacity = 1.0.
 *
 * Finally we call pruneDegreeOne(...) (if you have that).
 */

/*





void Graph::print() const {
    for (int i = 0; i < m_iNumNodes; ++i) {
        std::cout << "Node " << i << ": ";
        for (const auto& edge : m_adj[i]) {
            std::cout << "(" << edge->target << ", " << edge->capacity << ") ";
        }
        std::cout << std::endl;
    }
}

void Graph::readGraph(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int n = 0, m = 0;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'p') {
            std::string tmp;
            iss >> tmp >> tmp >> n >> m;
            break;
        }
    }

    if (n == 0) {
        throw std::runtime_error("Invalid or missing problem line in DIMACS file");
    }

    this->m_iNumNodes = 0;
    this->m_adj.resize(n);
    m_vertices.resize(this->m_iNumNodes);
    std::iota(m_vertices.begin(), m_vertices.end(), 0); // fill vertices with 0, 1, ..., maxNodeIdSeen

    do {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'a') {
            char dummy;
            int u, v;
            double cap = 1.0;
            iss >> dummy >> u >> v;
            if (!(iss >> cap)) cap = 1.0;
            this->addEdge(u - 1, v - 1, cap);  // convert to 0-based
        }
    } while (std::getline(infile, line));


}

const DiGraph& Graph::GetDiGraph() const{
    if(!m_diGraph) {
        auto temp = std::make_unique<DiGraph>(m_iNumNodes);
        for (int i = 0; i < m_iNumNodes; ++i) {
            for(const auto& edge : m_adj[i]) {
                if(i < edge->target) {
                    // To avoid adding the same undirected arc twice
                    temp->addUndirectedArc(i, edge->target, edge->capacity);
                }
            }
        }
        m_diGraph = std::move(temp);
    }
    return *m_diGraph;

}


void Graph::createDistanceMatrix() {
    // initialize the border elements from the matrix with the edge distances
    m_distanceMatrix.resize(m_iNumNodes, std::vector<double>(m_iNumNodes, std::numeric_limits<double>::infinity()));
    m_ArcMatrix.resize(m_iNumNodes, std::unordered_map<int, std::shared_ptr<Edge>>());

    for(int i = 0; i < m_iNumNodes; i++) {
        for(auto & edge : m_adj[i]) {
            {
                m_distanceMatrix[i][edge->target] = edge->distance;
                m_ArcMatrix[i][edge->target] = (edge);
            }
        }
    }

    // set the diagonal elements to 0
    for(int i = 0; i < m_iNumNodes; i++) {
        m_distanceMatrix[i][i] = 0.0;
        m_ArcMatrix[i][i] = nullptr; // No arc to self
    }

    for(int k = 0; k<m_iNumNodes; k++) {
        for(int i = 0; i<m_iNumNodes; i++) {
            double d_ik = m_distanceMatrix[i][k];

            for(int j = 0; j<m_iNumNodes; j++) {
                double d_ij = m_distanceMatrix[i][j];
                double d_kj = m_distanceMatrix[k][j];
                if(d_ij > d_ik + d_kj) {
                    m_distanceMatrix[i][j] = d_ik + d_kj;
                    m_ArcMatrix[i][j] = m_ArcMatrix[i][k];
                }
            }
        }
    }

}

const std::vector<std::vector< double> >& Graph::GetDistanceMatrix() const {
    return m_distanceMatrix;
}
const std::unordered_map<int, std::shared_ptr<Edge>>& Graph::GetArcMatrix(int node) const {
    return m_ArcMatrix[node];
}

double Graph::GetGraphDiameter() const {
    if(m_distanceMatrix.empty()) {
        throw std::runtime_error("Distance matrix is not initialized. Call createDistanceMatrix() first.");
        return -1;
    }

    double diameter = 0.0;
    for(int i = 0; i < m_iNumNodes; i++) {
        for(const auto& distance : m_distanceMatrix[i]) {
            if(distance > diameter) {
                diameter = distance;
            }
        }
    }
    return diameter;
}


// TODO: shortest path is not returned... for the example from 2 to 5...
std::vector<std::shared_ptr<Edge>> Graph::getShortestPath(int source, int target) const {
    const double INF = std::numeric_limits<double>::infinity();

    // static, so it can be returned by reference
    std::vector<std::shared_ptr<Edge>> path;


    std::vector<double> dist(m_iNumNodes, INF);
    std::vector<int> prev(m_iNumNodes, -1);
    dist[source] = 0.0;

    using P = std::pair<double, int>;
    std::priority_queue<P, std::vector<P>, std::greater<>> pq;
    pq.emplace(0.0, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        for (const auto& edge : m_adj[u]) {
            int v = edge->target;
            double weight = edge->distance;

            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    if (dist[target] == INF) {
        // No path exists
        return path;
    }

    // Reconstruct path from source to target
    std::vector<int> reverse_node_path;
    for (int u = target; u != -1; u = prev[u])
        reverse_node_path.push_back(u);
    std::reverse(reverse_node_path.begin(), reverse_node_path.end());

    for (size_t i = 0; i + 1 < reverse_node_path.size(); ++i) {
        int u = reverse_node_path[i];
        int v = reverse_node_path[i + 1];
        for (const auto& edge : m_adj[u]) {
            if (edge->target == v) {
                path.push_back(edge);
                break;
            }
        }
    }

    return path;
}

////////////////////////////////////////
//////////////// DiGraph ///////////////
////////////////////////////////////////


int DiGraph::numNodes() const {
    return m_iNumNodes;
}

int DiGraph::numArcs() const {
    return m_iNumArcs;
}

const Arc& DiGraph::GetArcById(int id) const{
    const auto& arc = m_IdToArc.find(id);
    return arc->second;
}

void DiGraph::addUndirectedArc(int source, int target, double capacity) {
    if (source < 0 || source >= m_iNumNodes || target < 0 || target >= m_iNumNodes) {
        throw std::out_of_range("Invalid node index");
    }

    Arc arc;
    arc.source = source;
    arc.target = target;
    arc.capacity = capacity;
    arc.id = this->m_idArc++;

    Arc rev_arc;
    rev_arc.source = target;
    rev_arc.target = source;
    rev_arc.capacity = capacity;
    rev_arc.id = this->m_idArc++;

    arc.rev_id = rev_arc.id;
    rev_arc.rev_id = arc.id;

    m_IdToArc[arc.GetId()] = arc;
    m_IdToArc[rev_arc.GetId()] = rev_arc;

    this->m_iNumArcs += 2; // Each undirected arc adds two directed arcs (one for each direction)
}


bool DiGraph::isValid() {
    bool result = true;
    for(const auto& [id, e] : GetArcs()) {
        if(id == e.GetId() ) {
            int rev_id = e.GetReverseArcId();
            Arc revArc = m_IdToArc[rev_id];
            if(e.GetId() < rev_id && e.GetId() + 1 != rev_id) {
                std::cout << "Invalid Edge/Rev Id" << std::endl;
                std::cerr << e.GetId() << " " <<  rev_id << std::endl;
                result &= false;
            }

            if(id % 2 == 0) {
                if(id + 1 != rev_id) {
                    std::cout << "Wrong enumeration" << std::endl;
                    std::cerr << e.GetId() << " " <<  rev_id << std::endl;
                    result &= false;
                }
            }else{
                if(id -1 != rev_id) {
                    std::cout << "Wrong enumeration" << std::endl;
                    std::cerr << e.GetId() << " " <<  rev_id << std::endl;
                    result &= false;
                }
            }


            if(revArc.GetId() != e.GetReverseArcId() || revArc.GetReverseArcId() != e.GetId()) {
                std::cout << "Error: " << std::endl;
                result &= false;
            }
        }else{
            result &= false;
        }
    }


    return result;
}

bool DiGraph::HasArcId(int id) const {

    if(m_IdToArc.find(id) != m_IdToArc.end()) {
        return true;
    }else{
        return false;
    }
}

void DiGraph::print() const {
    for(const auto& [id, arc] : this->GetArcs()) {
        std::cout << "Arc ID: " << id << "/" << arc.GetId() << " revID: " << arc.GetReverseArcId() << " ("  << arc.GetSource() << " , " << arc.GetTarget() << ", " << arc.GetCapacity() << ") " << std::endl;
    }
}

const std::unordered_map<int, Arc>& DiGraph::GetArcs() const {
    return this->m_IdToArc;
}

void DiGraph::readDiGraph(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int n = 0, m = 0;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'p') {
            std::string tmp;
            iss >> tmp >> tmp >> n >> m;
            break;
        }
    }

    if (n == 0) {
        throw std::runtime_error("Invalid or missing problem line in DIMACS file");
    }

    this->m_iNumNodes = n;
    // TODO: add here the functionality to fill m_vertices with 0, 1, ..., maxNodeIdSeen for GetVertices()

    do {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'a') {
            char dummy;
            int u, v;
            double cap = 1.0;
            iss >> dummy >> u >> v;
            if (!(iss >> cap)) cap = 1.0;
            this->addUndirectedArc(u - 1, v - 1, cap);  // convert to 0-based

        }
    } while (std::getline(infile, line));


}

 */
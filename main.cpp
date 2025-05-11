#include <iostream>
#include "src/graph.h"
#include "test/randomgraphs/random_graph_generator.h"
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    Graph G = RandomGraphGenerator::generate(10, 20, 1.0, 10.0);

    G.print();
    std::cout << "Number of nodes: " << G.numNodes() << std::endl;
    return 0;
}

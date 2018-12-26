/*
 * Name: Jerry Qian
 * PID: A14457025
 */

import java.util.Collection;
import java.util.Hashtable;
import java.util.List;
import java.util.PriorityQueue;
import java.util.*;

/**
 * Graph class which includes 4 different traversal algorithms.
 */

public class Graph {
    //storing graph as hash table with name as key and vertex as value
    Hashtable<String, Vertex> graph;

    /**
     * Constructor for Graph
     */
    public Graph() {
        graph = new Hashtable<>();
    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if given vertex
     * already exist in the graph.
     *
     * @param v vertex to be added to the graph
     * @throws IllegalArgumentException if two vertices with the same name are added.
     */
    public void addVertex(Vertex v) throws IllegalArgumentException {
        if (graph.containsKey(v.name)) {
            throw new IllegalArgumentException();
        }
        graph.put(v.name, v);
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {
        return graph.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        return graph.get(name);
    }

    /**
     * Adds a directed edge from vertex u to vertex v, Throws IllegalArgumentException if one of
     * the vertex does not exist
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight weight of the edge between vertex u and v
     * @throws IllegalArgumentException if one of the vertex does not exist
     */
    public void addEdge(String nameU, String nameV, Double weight) throws IllegalArgumentException {
        if (getVertex(nameU) == null || getVertex(nameV) == null) {
            throw new IllegalArgumentException();
        }

        Vertex u = getVertex(nameU);
        Vertex v = getVertex(nameV);
        //instantiate an edge
        Edge edge = new Edge(u, v, weight);
        //add edge but only to starting vertex
        u.edges.add(edge);
    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight  weight of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double weight) throws IllegalArgumentException{
        if (getVertex(nameU) == null || getVertex(nameV) == null) {
            throw new IllegalArgumentException();
        }

        Vertex u = getVertex(nameU);
        Vertex v = getVertex(nameV);
        //instantiate two edges
        Edge UV = new Edge(u, v, weight);
        Edge VU = new Edge(v, u, weight);
        //adds edges to both vertices
        u.edges.add(UV);
        v.edges.add(VU);
    }

    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
        return Math.sqrt(((vy - uy) * (vy - uy)) + ((vx - ux) * (vx - ux)));
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {
        for (Vertex uVertex: graph.values()) {
            for (Edge edge: uVertex.edges) {
                //obtain source and target vertices
                Vertex source = edge.source;
                Vertex target = edge.target;
                edge.distance = computeEuclideanDistance(source.x, source.y, target.x, target.y);
            }
        }
    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {
        for (Vertex vertex: graph.values()) {
            vertex.visited = false;
            //reset edges
            vertex.previous = null;
        }
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using DFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void DFS(String s, String t) {
        resetAllVertices();

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        //use stack to hold vertices
        Stack<Vertex> vertices = new Stack<>();
        //add first vertex to stack
        vertices.push(start);


        while (!vertices.isEmpty()) {
            //retrieve first vertex
            Vertex current = vertices.pop();

            //when all vertices up to the end have been visited, exit
            if (current.equals(end)) {
                return;
            }

            //if vertex has not been visited
            if (!current.visited) {
                //visit the vertex
                current.visited = true;
                for (Edge edge: current.edges) {
                    //find adjacent vertices
                    Vertex adjacent = edge.target;
                    if (!adjacent.visited) {
                        //set the previous edge
                        adjacent.previous = edge;
                        vertices.push(adjacent);
                    }
                }
            }
        }

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using BFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void BFS(String s, String t) {
        resetAllVertices();

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        //use linked list as queue
        LinkedList<Vertex> vertices = new LinkedList<>();

        //add first vertex to linked list (queue)
        vertices.add(start);

        while (!vertices.isEmpty()) {
            //removing the first element
            Vertex current = vertices.remove();

            //when all vertices up to the end have been visited, exit
            if (current.equals(end)) {
                return;
            }

            for (Edge edge: current.edges) {
                //obtain adjacent vertices
                Vertex adjacent = edge.target;
                //if vertex has not been visited
                if (!adjacent.visited) {
                    //set the previous edge
                    adjacent.previous = edge;
                    vertices.add(adjacent);
                    //visit the vertex
                    adjacent.visited = true;
                }
            }
        }

    }

    /**
     * Helper class for Dijkstra and A*, used in priority queue
     */
    private class CostVertex implements Comparable<CostVertex> {
        double cost;
        Vertex vertex;

        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
        }

        public int compareTo(CostVertex o) {
            return Double.compare(cost, o.cost);
        }
    }

    /**
     * Find the shortest path from vertex with name s to vertex with name t, using Dijkstra
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void Dijkstra(String s, String t) {
        resetAllVertices();

        //make a priority queue of cost vertices
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        //maps each vertex to its corresponding cost vertex in hashtable
        Hashtable<Vertex, CostVertex> costGraph = new Hashtable<>();

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        //creating a new graph consisting of cost vertices
        for (Vertex vertex : graph.values()) {
            CostVertex costVertex = new CostVertex(Double.POSITIVE_INFINITY, vertex);
            costGraph.put(vertex, costVertex);
        }
        //set the starting cost vertex's distance to 0
        costGraph.get(start).cost = 0;

        //offer all cost vertices into priority queue
        for (CostVertex costVertex: costGraph.values()) {
            queue.offer(costVertex);
        }

        //remove first cost vertex
        Vertex current = queue.poll().vertex;

        //loop through all vertices in queue
        while (!queue.isEmpty()) {
            //visit the first vertex
            current.visited = true;

            //exit the loop when end has been found
            if (current.equals(end)) {
                return;
            }

            //searching through adjacent vertices
            for (Edge e : current.edges) {
                //retrieve adjacent vertex and adjacent cost vertex
                Vertex adjacent = e.target;
                CostVertex adjacentCost = costGraph.get(adjacent);

                //calculating adjusted distance
                double adjustedCost = costGraph.get(current).cost + e.distance;

                //if vertex has not been visited and current cost is greater than adjusted cost
                //then update cost
                if (!adjacent.visited && adjacentCost.cost > adjustedCost) {
                    adjacentCost.cost = adjustedCost;
                    //keeping track of the edge it came from
                    adjacent.previous = e;

                    //remove vertex from queue and then reinsert it to maintain order
                    queue.remove(adjacentCost);
                    queue.offer(adjacentCost);
                }
            }
            //looping through next vertex
            current = queue.poll().vertex;
        }

    }

    /**
     * Helper method to calculate the h value in A*
     *
     * @param cur the current vertex being explored
     * @param goal the goal vertex to reach
     * @return the h value of cur and goal vertices
     */
    private double hValue(String cur, String goal) {
        Vertex current = getVertex(cur);
        Vertex end = getVertex(goal);
        return computeEuclideanDistance(end.x, end.y, current.x, current.y);

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void AStar(String s, String t) {
        resetAllVertices();

        //make a priority queue of cost vertices
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        //maps each vertex to its corresponding cost vertex in hashtable
        Hashtable<Vertex, CostVertex> costGraph = new Hashtable<>();

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        //creating a new graph consisting of cost vertices
        for (Vertex vertex : graph.values()) {
            CostVertex costVertex = new CostVertex(Double.POSITIVE_INFINITY, vertex);
            costGraph.put(vertex, costVertex);
        }
        //set the starting cost vertex's distance to 0
        costGraph.get(start).cost = 0;

        //offer all cost vertices into priority queue
        for (CostVertex costVertex: costGraph.values()) {
            queue.offer(costVertex);
        }

        //remove first cost vertex
        Vertex current = queue.poll().vertex;

        //loop through all vertices in queue
        while (!queue.isEmpty()) {
            //visit the first vertex
            current.visited = true;

            //exit the loop when end has been found
            if (current.equals(end)) {
                return;
            }

            //searching through adjacent vertices
            for (Edge e : current.edges) {
                //retrieve adjacent vertex and adjacent cost vertex
                Vertex adjacent = e.target;
                CostVertex adjacentCost = costGraph.get(adjacent);

                //calculating adjusted distance
                double adjustedCost = costGraph.get(current).cost + e.distance;

                //if vertex has not been visited and current cost is greater than adjusted cost
                //then update cost
                if (!adjacent.visited && adjacentCost.cost > adjustedCost) {
                    adjacentCost.cost = adjustedCost;
                    //keeping track of the edge it came from
                    adjacent.previous = e;

                    //remove vertex from queue and then reinsert it to maintain order
                    queue.remove(adjacentCost);
                    queue.offer(new CostVertex(adjacentCost.cost + hValue(adjacent.name, t), adjacent));
                }
            }
            //looping through next vertex
            current = queue.poll().vertex;
        }
    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) {

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);
        //create path to hold list of edges
        ArrayList<Edge> path = new ArrayList<>();

        Vertex current = end;
        //traverse from end vertex to start vertex
        while (!start.equals(current)) {
            //obtain previous edge
            Edge prev = current.previous;
            path.add(prev);
            //obtain previous vertex
            current = current.previous.source;
        }
        //reveres the list of edges
        Collections.reverse(path);
        return path;
    }

}

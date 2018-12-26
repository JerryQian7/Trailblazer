/*
 * Name: Jerry Qian
 * PID: A14457025
 */

import java.util.ArrayList;
import java.util.Hashtable;

/**
 * Designates a vertex within the graph.
 */
public class Vertex {

    public String name; // the name of this vertex
    public int x; // the x coordinates of this vertex on map
    public int y; // the y coordinates of this vertex on map
    public ArrayList<Edge> edges; //list of edges that this vertex has
    public boolean visited; //whether this vertex has been visited or not
    public Edge previous; //stores the previous edge it came from

    public Vertex(String name, int x, int y) {
        this.name = name;
        this.x = x;
        this.y = y;
        edges = new ArrayList<>();
        visited = false;
    }

    /**
     * Hashing function for the name of the vertex
     * @return the hashcode for the name
     */
    @Override
    public int hashCode() {
        // we assume that each vertex has a unique name
        return name.hashCode();
    }

    /**
     * Method that checks whether vertices are equal to each other
     * @param o vertex
     * @return whether they are equal
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (!(o instanceof Vertex)) {
            return false;
        }
        Vertex oVertex = (Vertex) o;

        return name.equals(oVertex.name) && x == oVertex.x && y == oVertex.y;
    }

    /**
     * Prints out vertex with its coordinates
     * @return vertex with its coordinates
     */
    public String toString() {
        return name + " (" + x + ", " + y + ")";
    }

}
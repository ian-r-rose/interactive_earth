/************************ /
    regular_grid.h
Define an abstract interface to a simple cartesian mesh,
with ways to iterate over the mesh cells, as well as query
for indexing and geometric information.  For use with
cartesian grids in serial finite difference
calculations
/ ***********************/

#include <iterator>
#include <cmath>

#ifndef REGULAR_GRID_H
#define REGULAR_GRID_H


//Struct for holding a 2D point.
//Can refer both to x,y and theta,r,
//though the memory layout is the same.
struct Point
{
  union
  {
    struct{ double x,y; };
    struct{ double theta, r;};
  };
};

inline int fast_floor(double x)
{
  int i = (int)x;  //truncate
  return i - ( i > x );  //handle negatives
}

class RegularGrid
{
  public:

    class iterator;
    class reverse_iterator;

    //Basic cell data type, which knows its own index,
    //geometric properties, and how to access its neighbors.
    class Cell
    {
      public:

        Cell (unsigned int i, const RegularGrid &g): id(i), grid(g) {}
        Cell (Cell &c): id(c.id), grid(c.grid) {}
        Cell &operator=(const Cell &rhs) { id=rhs.id; return *this;}

        //get the x,y index of the cell
        int xindex() { return id % grid.nx; }
        int yindex() { return id / grid.nx; }
        //get the r,theta index of the cell: identical to xindex, yindex
        int thetaindex() { return id % grid.nx; }
        int rindex() { return id / grid.nx; }

        //get grid indices of the neighbors
        int left() { return (id % grid.nx == 0 ? id+grid.nx-1 : id-1); }
        int right() { return ( (id+1)%grid.nx == 0 ? id-grid.nx+1 : id+1); }
        int up() { return (id + grid.nx >= grid.ncells ? id-grid.nx*(grid.ny-1) : id + grid.nx); }
        int down() { return ( id-grid.nx < 0 ? id+grid.nx*(grid.ny-1) : id-grid.nx) ; }
        int upleft() { return id + (id + grid.nx >= grid.ncells ? -grid.nx*(grid.ny-1) : grid.nx) + (id % grid.nx == 0 ? grid.nx-1 : -1); }
        int upright() { return id + (id + grid.nx >= grid.ncells ? -grid.nx*(grid.ny-1) : grid.nx) + ( (id+1)%grid.nx == 0 ? -grid.nx+1 : 1); }
        int downleft() { return id + (id-grid.nx < 0 ? grid.nx*(grid.ny-1) : -grid.nx) + ( id%grid.nx == 0 ? grid.nx-1 : -1); }
        int downright() { return id + (id+grid.nx <0 ? grid.nx*(grid.ny-1) : -grid.nx) + ( (id+1)%grid.nx == 0 ? -grid.nx+1 : 1); }
        int self() { return id; };

        //query for boundary information
        bool at_top_boundary() {return (id + grid.nx >= grid.ncells);}
        bool at_bottom_boundary() {return (id-grid.nx < 0);}
        bool at_left_boundary() {return (id%grid.nx == 0);}
        bool at_right_boundary() {return (id+1)%grid.nx == 0;}
        bool at_boundary() {return at_top_boundary() || at_bottom_boundary() || at_left_boundary() || at_right_boundary(); }

        //Get some location information.  This is more complicated than a nonstaggered grid,
        //as different of the properties will be found on different parts of the cell.
        Point location() { Point p; p.x = xindex()*grid.dx; p.y = yindex()*grid.dy;  return p;} //lower left
        double radius() { return grid.r_inner + (rindex()*grid.dr); }


      private:
        const RegularGrid& grid;  //Const reference to the grid
        int id; //id of the cell, which will correspond to its index in vectors

      friend class RegularGrid::iterator;
      friend class RegularGrid::reverse_iterator;

    };

    //Lightweight iterator that starts at the first grid cell
    //(bottom left), and iterates over all the cells.
    class iterator : public std::iterator<std::input_iterator_tag, Cell>
    {
      private:
        int id;
        const RegularGrid& grid;
        Cell c;
      public:
        iterator(int i, const RegularGrid& g): id(i), grid(g), c(id,grid) {};
        iterator(const iterator &it) : id(it.id), grid(it.grid), c(id, grid) {};
        iterator& operator++() { ++id; ++c.id; return (*this);}
        iterator& operator++(int) { return operator++();}
        bool operator==(const iterator &rhs) {return id == rhs.id; };
        bool operator!=(const iterator &rhs) {return id != rhs.id; };
        Cell &operator*() {return c;}
        Cell* operator->() {return &c;}
    };
    //Lightweight iterator that starts at the last grid cell (top right),
    //and iterates backwards over all the cells.
    class reverse_iterator : public std::iterator<std::input_iterator_tag, Cell>
    {
      private:
        int id;
        const RegularGrid& grid;
        Cell c;
      public:
        reverse_iterator(int i, const RegularGrid& g): id(i), grid(g), c(id,grid) {};
        reverse_iterator(const reverse_iterator &it) : id(it.id), grid(it.grid), c(id, grid) {};
        reverse_iterator& operator++() { --id; --c.id; return (*this);}
        reverse_iterator& operator++(int) { return operator++();}
        bool operator==(const reverse_iterator &rhs) {return id == rhs.id; };
        bool operator!=(const reverse_iterator &rhs) {return id != rhs.id; };
        Cell &operator*() {return c;}
        Cell* operator->() {return &c;}
    };

    //const members so we can query directly
    const union { double lx; double ltheta; }; //length in x/theta direction
    const union { double ly; double lr; }; //length in y/r direction
    const union { int nx; int ntheta; }; //Numer of nodes in x/theta direction
    const union { int ny; int nr; }; //Numer of nodes in y/r direction
    const union { double dx; double dtheta; }; //node spacing in x/theta direction
    const union { double dy; double dr; }; //node spacing in y/r direction
    const int ncells; //Total number of cells
    const double r_inner, r_outer;

    RegularGrid(const double inner_radius, const unsigned int num_theta, const unsigned int num_r)
                  : r_inner(inner_radius), r_outer(1.0), ltheta(2.0*M_PI), lr(1.0-inner_radius),
                    ntheta(num_theta), nr(num_r), dtheta(ltheta/ntheta), dr(lr/(nr-1)), ncells(ntheta*nr)
                  {}
    const iterator begin() { return iterator(0, *this); }
    const iterator end() {return iterator(ncells, *this);}
    const reverse_iterator rbegin() { return reverse_iterator(ncells-1, *this); }
    const reverse_iterator rend() {return reverse_iterator(-1, *this);}

    //Get handles to cell id and cell iterators at a point
    inline int cell_id( const Point &p) const { int xindex = fast_floor(p.x/dx); int yindex=fast_floor(p.y/dy);
                                   return keep_in_domain(xindex, yindex); }
    inline int keep_in_domain( int xindex, int yindex) const { return nx*(yindex < 0 ? 0 : (yindex >= ny ? ny-1: yindex))
                              + (xindex % nx + nx) %nx ; }
    iterator cell_at_point( const Point &p) { return iterator(cell_id(p), *this); }

    iterator lower_left_cell( const Point &p) { return cell_at_point(p); }

};


//Cubic lagrange interpolation at an arbitrary point, given the values of the
//function at a nine-point stencil around it.  I am currently using linear
//rather than cubic interpolation for performance reasons.
inline double cubic_interp_2d( double x, double y, double ul, double u, double ur,
                                  double l, double c, double r, double dl, double d, double dr)
{
  return   ul*(x)*(x-1.0)*(y)*(y+1.0)/4.0
         - u*(x-1.0)*(x+1.0)*(y)*(y+1.0)/2.0
         + ur*(x+1.0)*(x)*(y)*(y+1.0)/4.0
         - l*(x)*(x-1.0)*(y-1.0)*(y+1.0)/2.0
         + c*(x-1.0)*(x+1.0)*(y-1.0)*(y+1.0)
         - r*(x+1.0)*(x)*(y-1.0)*(y+1.0)/2.0
         + dl*(x)*(x-1.0)*(y)*(y-1.0)/4.0
         - d*(x-1.0)*(x+1.0)*(y)*(y-1.0)/2.0
         + dr*(x)*(x+1.0)*(y)*(y-1.0)/4.0;
}

//Linear interpolation at a point, given the values at a four-point stencil
//around it.
inline double linear_interp_2d(double x, double y, double ul, double ur, double dl, double dr)
{
  return - ul * (x-1.0) * (y)
         + ur * (x) * (y)
         + dl * (x-1.0) * (y-1.0)
         - dr * (x) * (y-1.0);
}

#endif

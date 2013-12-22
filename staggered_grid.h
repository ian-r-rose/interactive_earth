/************************ /
    StaggeredGrid.h
Define an abstract interface to a simple cartesian mesh,
with ways to iterate over the mesh cells, as well as query
for indexing and geometric information.  For use with
staggered cartesian grids in serial finite difference 
calculations
/ ***********************/

#include <iterator>
#include <iostream>

struct Point
{
  double x;
  double y;
};

class StaggeredGrid
{
  public:

    class iterator;

    class Cell 
    {
      public:

        Cell (unsigned int i, const StaggeredGrid &g): id(i), grid(g) {};
        Cell (Cell &c): id(c.id), grid(c.grid) {};
        Cell &operator=(const Cell &rhs) { id=rhs.id; return *this;}
 
        //get the x,y index of the cell
        int xindex() { return id % grid.nx; };
        int yindex() { return id / grid.nx; };
        
        //get grid indices of the neighbors
        int left() { return (id % grid.nx == 0 ? -1 : id-1); };
        int right() { return ( (id+1)%grid.nx == 0 ? -1 : id+1); };
        int up() { return (id + grid.nx >= grid.ncells ? -1 : id + grid.nx); };
        int down() { return ( id-grid.nx < 0 ? -1 : id-grid.nx) ; };
        int self() { return id; };
 
        //query for boundary information
        bool at_top_boundary() {return (id + grid.nx >= grid.ncells);};
        bool at_bottom_boundary() {return (id-grid.nx < 0);};
        bool at_left_boundary() {return (id%grid.nx == 0);};
        bool at_right_boundary() {return (id+1)%grid.nx == 0;};
        bool at_boundary() {return at_top_boundary() || at_bottom_boundary() || at_left_boundary() || at_right_boundary(); }

        //Get some location information
        Point center() { Point p; p.x = xindex()*grid.dx + grid.dx/2.0; p.y = yindex()*grid.dy + grid.dy/2.0;  return p;};
        Point corner() { Point p; p.x = xindex()*grid.dx; p.y = yindex()*grid.dy;  return p;};
        Point hface() { Point p; p.x = xindex()*grid.dx + grid.dx/2.0; p.y = yindex()*grid.dy;  return p;};
        Point vface() { Point p; p.x = xindex()*grid.dx; p.y = yindex()*grid.dy + grid.dy/2.0;  return p;};
        
      private:
        const StaggeredGrid& grid;
        int id;

      friend class StaggeredGrid::iterator;

    };

    class iterator : public std::iterator<std::input_iterator_tag, Cell>
    {
      private:
        int id;
        const StaggeredGrid& grid;
        Cell c;
      public:
        iterator(unsigned int i, const StaggeredGrid& g): id(i), grid(g), c(id,grid) {};
        iterator(const iterator &it) : id(it.id), grid(it.grid), c(id, grid) {};
        iterator& operator++() { ++id; ++c.id; return (*this);}
        iterator& operator++(int) { return operator++();}
        bool operator==(const iterator &rhs) {return id == rhs.id; };
        bool operator!=(const iterator &rhs) {return id != rhs.id; };
        Cell &operator*() {return c;}
        Cell* operator->() {return &c;}
    };
        
    //const members so we can query directly
    const double lx, ly; //Length in x,y directions
    const int nx, ny; //Number of cells in x,y directions
    const double dx, dy; //Number of cells in x,y directions
    const int ncells;

    StaggeredGrid(const double lenx, const double leny, const unsigned int numx, const unsigned int numy)
                  : lx(lenx), ly(leny), nx(numx), ny(numy), dx(lx/nx), dy(ly/ny), ncells(nx*ny) {}; 
    const iterator begin() { return iterator(0, *this); };
    const iterator end() {return iterator(ncells, *this);};
};

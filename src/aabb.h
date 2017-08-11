#ifndef AABB_H
#define AABB_H

#include <algorithm>
#include <stack>
#include <cassert>
#include <vector>

using namespace std;

namespace path_planning
{

struct AABBNode
{
  // aabb fields
  double max_x_;
  double min_x_;
  double max_y_;
  double min_y_;
  // time dimension
  int min_t_;
  int max_t_;

  // node fields
  int parent_ = -1;
  int left_   = -1;
  int right_  = -1;

  AABBNode(
     double  min_x
    , double max_x
    , double  min_y
    , double  max_y
    , int  min_t
    , int  max_t
    ){
    assert(min_x < max_x);
    assert(min_y < max_y);
    assert(min_t <= max_t);
    min_x_ = min_x;
    max_x_ = max_x;
    min_y_ = min_y;
    max_y_ = max_y;
    min_t_ = min_t;
    max_t_ = max_t;
  }

  static AABBNode FromPoint
    (double x
    , double y
    , int t
    , double width
    , double height){
  
    double half_width = 0.5 * width;
    double half_height = 0.5 * height;
    return AABBNode(x - half_width
                  , x + half_width
                  , y - half_height
                  , y + half_height
                  , t, t);
  }

  double Width(){
    return max_x_ - min_x_;
  }

  double Height(){
    return max_y_ - min_y_;
  }

  vector<double> Center(){
    return {min_x_ + 0.5 * Width(), min_y_ + 0.5 * Height()};
  }

  bool IsOverlap(const AABBNode& other){
    return 
         min_x_ <= other.max_x_ 
      && max_x_ >= other.min_x_ 
      && min_y_ <= other.max_y_
      && max_y_ >= other.min_y_
      && min_t_ <= other.max_t_
      && max_t_ >= other.min_t_;
  }

  const friend AABBNode operator + (const AABBNode& left, const AABBNode& right){
    return AABBNode(
              min(left.min_x_, right.min_x_)
            , max(left.max_x_, right.max_x_)
            , min(left.min_y_, right.min_y_)
            , max(left.max_y_, right.max_y_)
            , min(left.min_t_, right.min_t_)
            , max(left.max_t_, right.max_t_)
            );
  }

  friend std::ostream& operator << (std::ostream& os, const AABBNode& node) {
      os << "AABBNode=["
      << "min_x: " << node.min_x_  << ", "
      << "max_x: " << node.max_x_  << ", "
      << "min_y: " << node.min_y_  << ", "
      << "max_y: " << node.max_y_  << ", "
      << "min_t: " << node.min_t_  << ", "
      << "max_t: " << node.max_t_  << " | "
      << "parent: "<< node.parent_ << ", "
      << "left: "  << node.left_   << ", "
      << "right: " << node.right_  << "]";
      return os;
  }

  bool IsLeaf() { 
    return left_ == -1 && right_ == -1; 
  }
};


struct AABBTree {
private:
  enum node_side {LEFT, RIGHT};
public:
  vector<AABBNode> nodes_;
  int root_ = -1;

  static AABBTree Build(vector<double> xs, vector<double> ys, double width, double height){
    assert(xs.size() > 0);
    assert(xs.size() == ys.size());

    AABBTree tree;
    double half_width = width/2;
    double half_height = height/2;

    stack<vector<int>> tree_stack;
    tree_stack.push({-1, 0, int(xs.size()), LEFT});

    while (!tree_stack.empty()){
      auto segment = tree_stack.top();

      int parent = segment[0];
      int start = segment[1];
      int end = segment[2];
      int side = segment[3];

      double min_x = *(min_element(xs.begin()+start, xs.begin()+end)) - half_width ;
      double max_x = *(max_element(xs.begin()+start, xs.begin()+end)) + half_width ;
      double min_y = *(min_element(ys.begin()+start, ys.begin()+end)) - half_height;
      double max_y = *(max_element(ys.begin()+start, ys.begin()+end)) + half_height;

      AABBNode new_node(min_x, max_x, min_y, max_y, start, end);
      new_node.parent_ = parent;

      tree.nodes_.push_back(new_node);
      
      int new_index = tree.nodes_.size()-1;
      if (parent != -1){
        if (side == LEFT){
          tree.nodes_[parent].left_ = new_index;  
        }      
        else{
          tree.nodes_[parent].right_ = new_index; 
        }       
      }
      else {
        tree.root_ = new_index;
      }

      tree_stack.pop();
      
      int diff = (end-start)/2;
      if (diff>0){
        int middle = start+diff;
        tree_stack.push({new_index, middle, end, RIGHT});
        tree_stack.push({new_index, start, middle, LEFT});
      }
    }

    return tree;
  }

  int CollisionIndex(AABBNode node){
    int result = -1;
    stack<int> tree_stack;
    tree_stack.push(root_);
    while (!tree_stack.empty()){
      int index = tree_stack.top();
      tree_stack.pop();
      if (index != -1){
        auto check_node = nodes_[index];
        if (node.IsOverlap(check_node)){
          if (check_node.IsLeaf()){
            result = index;
            break;
          }
          else{
            if (check_node.right_ != -1){
              tree_stack.push(check_node.right_);
            }
            if (check_node.left_ != -1){
              tree_stack.push(check_node.left_);
            }
          }
        }
      }
    }
    return result;
  }

  int CollisionIndex(vector<double> xs, vector<double> ys, double width, double height){
    assert(xs.size()==ys.size());
    int result = -1;
    for (int i=0; i<xs.size(); i++){
      auto node = AABBNode::FromPoint(xs[i], ys[i], i, width, height);
      int index = CollisionIndex(node);
      if (index != -1){
        result = index;
        break;
      }
    }
    return result;
  }

  bool IsCollideWith(AABBNode node){
    int index = CollisionIndex(node);
    return index != -1;
  }
  
  bool IsCollideWith(double x, double y, int t, double width, double height){
    auto node = AABBNode::FromPoint(x, y, t, width, height);
    return IsCollideWith(node);
  }

  bool IsCollideWith(vector<double> xs, vector<double> ys, double width, double height){
    int index = CollisionIndex(xs, ys, width, height);
    return index != -1;
  }


  friend std::ostream& operator << (std::ostream& os, const AABBTree& tree) {

    os << "AABBtree=[size: "<< tree.nodes_.size()<<", nodes: \n";

    if (tree.root_ != -1){
      stack<vector<int>> tree_stack;
      tree_stack.push({tree.root_, 0});
      while (!tree_stack.empty()){
        int index = tree_stack.top()[0];
        int level = tree_stack.top()[1];
        auto node = tree.nodes_[index];

        tree_stack.pop();
        for (int i=0; i<level; i++){
          cout<<"  ";
        }
        cout<<index<<":"<< node<<"\n";
        
        if (node.right_ != -1){
          tree_stack.push({node.right_, level + 1});        
        } 
        if (node.left_ != -1) {
          tree_stack.push({node.left_, level + 1});
        }
      }
    }
    os <<"]";
    return os;
  }

};

} /* path_planning */

#endif /* AABB_H */
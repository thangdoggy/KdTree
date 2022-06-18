#include <iostream>
#include <math.h>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>

// By QT 1/2022
// Implementation base on Chapter 9 - Advanced Data Structure and Algorithm - Marcello La Rocca

using namespace std;

auto start = chrono::steady_clock::now();

int LEVEL = 0;

int k;

class KdTree
{
public:
    class KdNode;
    KdNode *root;

    KdTree(int dimension)
    {
        k = dimension;
        root = nullptr;
    }
    // Taking an optional array of points as input
    KdTree(float **points, int numOfPoints, int dimension)
    {
        k = dimension;
        root = constructKdTree(points, 0, numOfPoints - 1, 0);
    }

    // Methods of K-d trees
    KdNode *constructKdTree(float **points, int left, int right, int level = 0);
    static bool mySort(float point1[], float point2[]);
    int partition(float **points, int left, int right, int level);

    KdNode *search(float target[]);
    KdNode *searchHelper(KdNode *node, float target[]);

    KdNode *insert(float newPoint[]);
    KdNode *insertHelper(KdNode *node, float newPoint[], int level = 0);

    KdNode *remove(float point[]);
    KdNode *removeHelper(KdNode *node, float point[]);
    KdNode *findMin(KdNode *node, int coordinateIndex);

    KdNode *nearestNeighbor(float target[]);
    KdNode *nearestNeighbor(KdNode *node, float target[], float &nnDist, KdNode *nn);
    float distance(float point1[], float point2[]);

    // Helper functions
    float getNodeKey(KdNode *node);
    float getPointKey(float point[], int level);
    int compare(float point[], KdNode *node);
    float splitDistance(float point[], KdNode *node);
    bool arePointsSame(float point1[], float point2[]);

    // print the tree structure for local testing
    void printBinaryTree(string prefix, const KdNode *root, bool isLeft, bool hasRightSibling)
    {
        if (!root && isLeft && hasRightSibling)
        {
            cout << prefix << "|--\n";
        }
        if (!root)
            return;
        cout << prefix;
        if (isLeft && hasRightSibling)
            cout << "|--";
        else
            cout << "|--";
        cout << "(";
        for (int i = 0; i < k - 1; i++)
        {
            cout << root->point[i] << ", ";
        }
        cout << root->point[k - 1] << ")" << endl;
        printBinaryTree(prefix + (isLeft && hasRightSibling ? "|  " : "   "), root->left, true, root->right);
        printBinaryTree(prefix + (isLeft && hasRightSibling ? "|  " : "   "), root->right, false, root->right);
    }

    void print()
    {
        printBinaryTree("", root, false, false);
    }

    class KdNode
    {
    public:
        float *point;
        KdNode *left;
        KdNode *right;
        int level;

        KdNode(float point[], KdNode *left, KdNode *right, int level)
        {
            this->point = new float[k];
            for (int i = 0; i < k; i++)
            {
                this->point[i] = point[i];
            }
            this->left = left;
            this->right = right;
            this->level = level;
        }
    };
};

// Methods of K-d trees //////////////////////////////////////////////////////////////

// Constructs a k-d tree from a set of points
KdTree::KdNode *KdTree::constructKdTree(float **points, int left, int right, int level)
{
    int size = right - left + 1;
    if (size == 0)
        return nullptr;
    else if (size == 1)
        return new KdNode(points[left], nullptr, nullptr, level);
    else
    {
        int median = left + partition(points, left, right, level);
        KdNode *leftTree = constructKdTree(points, left, median - 1, level + 1);
        KdNode *rightTree = constructKdTree(points, median + 1, right, level + 1);
        return new KdNode(points[median], leftTree, rightTree, level);
    }
}

bool KdTree::mySort(float point1[], float point2[])
{
    return point1[LEVEL % k] < point2[LEVEL % k];
}

int KdTree::partition(float **points, int left, int right, int level)
{
    LEVEL = level;
    sort(points + left, points + right, mySort);
    return (right - left) / 2;
}

/** Search returns the tree node that contains a target point,
 * if the point is stored in the tree; it returns null otherwise.
 */
KdTree::KdNode *KdTree::search(float target[])
{
    return searchHelper(root, target);
}

KdTree::KdNode *KdTree::searchHelper(KdNode *node, float target[])
{
    if (node == nullptr)
        return nullptr;
    else if (arePointsSame(node->point, target))
        return node;
    else if (compare(target, node) < 0)
        return searchHelper(node->left, target);
    else
        return searchHelper(node->right, target);
}

/** Inserts a point on the tree. The method will
 * return a pointer to the root of the (sub)tree
 * containing the new point.
 */
KdTree::KdNode *KdTree::insert(float newPoint[])
{
    return insertHelper(root, newPoint, 0);
}

KdTree::KdNode *KdTree::insertHelper(KdNode *node, float newPoint[], int level)
{
    if (node == nullptr)
        return new KdNode(newPoint, nullptr, nullptr, level);
    else if (arePointsSame(node->point, newPoint))
        return node;
    else if (compare(newPoint, node) < 0)
    {
        node->left = insertHelper(node->left, newPoint, node->level + 1);
        return node;
    }
    else
    {
        node->right = insertHelper(node->right, newPoint, node->level + 1);
        return node;
    }
}

/** Removes a point from the tree rooted at node.
 * The method will return the root of the tree after completing the operation.
 */
KdTree::KdNode *KdTree::remove(float point[])
{
    return removeHelper(root, point);
}

KdTree::KdNode *KdTree::removeHelper(KdNode *node, float point[])
{
    if (node == nullptr)
        return nullptr;
    else if (arePointsSame(node->point, point))
    {
        if (node->right != nullptr)
        {
            KdNode *minNode = findMin(node->right, node->level % k);
            KdNode *newRight = removeHelper(node->right, minNode->point);
            return new KdNode(minNode->point, node->left, newRight, node->level);
        }
        else if (node->left != nullptr)
        {
            KdNode *minNode = findMin(node->left, node->level % k);
            KdNode *newRight = removeHelper(node->left, minNode->point);
            return new KdNode(minNode->point, nullptr, newRight, node->level);
        }
        else
        {
            // delete node;
            return nullptr;
        }
    }
    else if (compare(point, node) < 0)
    {
        node->left = removeHelper(node->left, point);
        return node;
    }
    else
    {
        node->right = removeHelper(node->right, point);
        return node;
    }
}

// Finds the node in the tree with the minimum value for the coordinate at a given index
KdTree::KdNode *KdTree::findMin(KdNode *node, int coordinateIndex)
{
    if (node == nullptr)
        return nullptr;
    else if ((node->level % k) == coordinateIndex)
    {
        if (node->left == nullptr)
            return node;
        else
            return findMin(node->left, coordinateIndex);
    }
    else
    {
        KdNode *leftMin = findMin(node->left, coordinateIndex);
        KdNode *rightMin = findMin(node->right, coordinateIndex);
        if (leftMin == nullptr && rightMin == nullptr)
            return node;
        else if (leftMin == nullptr)
        {
            if (node->point[coordinateIndex] < rightMin->point[coordinateIndex])
                return node;
            else
                return rightMin;
        }
        else if (rightMin == nullptr)
        {
            if (node->point[coordinateIndex] < leftMin->point[coordinateIndex])
                return node;
            else
                return leftMin;
        }

        if (node->point[coordinateIndex] < leftMin->point[coordinateIndex] && node->point[coordinateIndex] < rightMin->point[coordinateIndex])
            return node;
        else if (leftMin->point[coordinateIndex] < rightMin->point[coordinateIndex])
            return leftMin;
        else
            return rightMin;
    }
}

/** Finds the closest point to a given target. We also pass
 * the best values found so far for nearest neighbor (NN)
 * and its distance to help pruning. These values default
 * to null, infinity for a call on the tree root: nnDist = inf, nn = null
 */
KdTree::KdNode *KdTree::nearestNeighbor(float target[])
{
    float nnDist = INT_MAX;
    KdNode *nn = nullptr;
    return nearestNeighbor(root, target, nnDist, nn);
}

KdTree::KdNode *KdTree::nearestNeighbor(KdNode *node, float target[], float &nnDist, KdNode *nn)
{
    if (node == nullptr)
        return nn;
    else
    {
        float dist = distance(node->point, target);
        if (dist < nnDist)
        {
            nnDist = dist;
            nn = node;
        }
        KdNode *closeBranch, *farBranch;
        if (compare(target, node) < 0)
        {
            closeBranch = node->left;
            farBranch = node->right;
        }
        else
        {
            closeBranch = node->right;
            farBranch = node->left;
        }
        nn = nearestNeighbor(closeBranch, target, nnDist, nn);
        if (splitDistance(target, node) < nnDist)
        {
            nn = nearestNeighbor(farBranch, target, nnDist, nn);
        }
        return nn;
    }
}

float KdTree::distance(float point1[], float point2[])
{
    float s = 0;
    for (int i = 0; i < k; i++)
    {
        s += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return sqrt(s);
}

// Helper functions

/** Given a tree node, returns the value of the
 * coordinate that needs to be used, given the
 * level at which the node is stored
 */
float KdTree::getNodeKey(KdNode *node)
{
    return getPointKey(node->point, node->level);
}

/** Given a point (a tuple with k values)
 * and an index for the level, returns the
 * tuple entry that should be used in
 * comparisons for nodes at that level
 */
float KdTree::getPointKey(float point[], int level)
{
    return point[level % k];
}

/** Compares a point to a node,
 * if the node’s point matches the first argument return -1 if even level and otherwise;
 * a value lower than 0 if the point is on the “left” of the node,
 * greater than 0 otherwise
 */
int KdTree::compare(float point[], KdNode *node)
{
    float s = getPointKey(point, node->level) - getNodeKey(node);
    if (s == 0)
        return node->level % 2 == 0 ? -1 : 1;
    else if (s > 0)
        return 1;
    else
        return -1;
}

/** Computes the distance between a point
 * and its projection on the split line passing through a node
 */
float KdTree::splitDistance(float point[], KdNode *node)
{
    return abs(getPointKey(point, node->level) - getNodeKey(node));
}

// Method to determine if two Points are same in K Dimensional space
bool KdTree::arePointsSame(float point1[], float point2[])
{
    // Compare individual point's values
    for (int i = 0; i < k; ++i)
        if (point1[i] != point2[i])
            return false;

    return true;
}

int main()
{
    // We are running nearest neighbor of (23.56,4.124,11.01)
    int numOfPoint;
    cout << "Number of point: ";
    cin >> numOfPoint;
    float **point3D = new float *[numOfPoint];
    for (int i = 0; i < numOfPoint; i++)
    {
        point3D[i] = new float[3];
    }

    random_device rd;
    mt19937 rng(rd());
    uniform_int_distribution<int> uni(0, 100000);

    for (int i = 0; i < numOfPoint; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            auto random_integer = uni(rng) / 1000.0;
            point3D[i][j] = random_integer;
        }
    }

    float target[3] = {23.56, 4.124, 11.01};
    // Construct a balance kd-tree
    KdTree *myTree = new KdTree(point3D, numOfPoint, 3);
    KdTree::KdNode *NN = myTree->nearestNeighbor(target);
    cout << "(" << NN->point[0] << ", " << NN->point[1] << ", " << NN->point[2] << ")\n";
    // myTree->print();

    cout << "\n";
    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << chrono::duration<double, milli>(diff).count() << " ms" << endl;
    return 0;
    /*
    KdTree* myTree = new KdTree(3);
    myTree->root = myTree->insert(point3D[0]);
    myTree->root = myTree->insert(point3D[1]);
    myTree->root = myTree->insert(point3D[3]);
    myTree->print(myTree->root);
    */
}

#include <iostream>
using namespace std;

void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

int main() {
    int myArray[5] = { 1, 2, 3, 4, 5 };
    printArray(myArray, 5);
    return 0;
}

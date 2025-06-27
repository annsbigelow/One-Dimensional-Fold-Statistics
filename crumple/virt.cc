#include <cstdio>

class animal {
    public:
        virtual void name() {
            puts("I am generic animal.");
        }
        virtual void sound() = 0; // Pure virtual function
};

class cat : public animal {
    public:
        virtual void name() {
            puts("I am a cat.");
        }
        virtual void sound() {
            puts("Meow");
        }
};

class polar_bear : public animal {
    public:
        virtual void name() {
            puts("I am a polar bear.");
        }
        virtual void sound() {
            puts("ROAAARRRRR!!!!");
        }
};

int main() {
    // This will give a compiler error because 'animal' has at
    // least one pure virtual function
    //animal a;

    polar_bear a;
    a.sound();
}

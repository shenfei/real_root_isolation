#include<iostream>
#include<string>

using namespace std;

int main()
{
    string line;
    do{
        getline(cin,line);
    }while (line.find("The order") == string::npos);

    while (getline(cin,line))
    {
        if (line.find('=') != string::npos)
            break;
        cout << line.substr(1) <<endl;
    }

    return 0;
}

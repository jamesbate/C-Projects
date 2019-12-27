#include <iostream>
#include <stdlib.h>
#include <ctime>
using namespace std;
int seed = 0;
//I use this as sometimes clocks ticks aren't fast enough to provide distinct seeds
unsigned int t = time(0);

bool one_sink_flag = false;
bool two_sink_flag = false;
bool victory_flag = false;
int argrid_size[2] = {10,10};
static int boat_number = 4;


class player
{
public:
    int ammo = 30;
    void receive_input();
    int strike_pos[2] = {0,0};
};

class boat
{
public:
    int arend_one[2];
    int arend_two[2];
    //delete
    int boat_pos[2][4] = {{0}};
    int boat_condition[4] = {{0}};  //1 indicates hit
    void boat_condition_set();
    short ndir;
    short nlength;
    bool afloat = true;
    void afloat_check();
};

void sink_message(boat boats[])
{
    short nsum = 0;
    for(short i=0;i<boat_number;i++)
    {
        if(boats[i].afloat == 0)
            nsum += 1;
    }
    if(nsum == 1 && one_sink_flag == 0)
    {
        cout << "One Boat Down!!!" << endl;
        one_sink_flag = 1;
    }
     if(nsum == 2 && two_sink_flag == 0)
    {
        cout << "Two Boats Down!!!" << endl;
        two_sink_flag = 1;
    }
    if(nsum == 3)
    {
        victory_flag = true;
        cout << "You WON!" << endl;
    }

}
void boat::boat_condition_set()
{
        for(int i = 0;i < boat_number;i++)
        {
            for(int j = 0;j<4-nlength;j++)
            boat_condition[4-j-1] = 1;
        }
    //configure condition array for boats according to their length
}


void boat::afloat_check()
{
    for(int i = 0;i<4;i++)
    {
        if(boat_condition[i] == 0)
        {
            return;
        }

    }
    afloat = false;
    return;
}

int rand_gen(int l, int u)
{
    //random number between l and u
    seed++;
    srand(seed + t);
    return rand()%(u-l+1)+l;
}

int rand_gen()
{
    seed++;
    srand(seed + t);
    return rand()%100+1;
    //overloaded to give default l, u vales
}


void player::receive_input()
{
    cout << "xpos:\t" << endl;
    cin >> strike_pos[0];
    cout << "ypos:\t" << endl;
    cin >> strike_pos[1];
    ammo--;
}
//get player decision

bool hit_result(player james,boat boats[])
{
    for(int k=0;k<boat_number;k++)//sum over boats
    {
            for(int j = 0;j<boats[k].nlength;j++)
            {
                if(james.strike_pos[0] == boats[k].boat_pos[0][j]&&james.strike_pos[1] == boats[k].boat_pos[1][j])
                   {
                       boats[k].boat_condition[j] = 1;
                       return true;
                   }
            }
    }
    return false;
}

int main()
{
    boat boats[boat_number];     //create three boat objects
    player james;

    for(short i = 0;i<boat_number;i++)
    {
        boats[i].arend_one[0] = rand_gen(1,argrid_size[0] - 2);
        boats[i].arend_one[1] = rand_gen(1,argrid_size[1] - 2);  //grid 0 to 9, placeable region 1,8
        boats[i].arend_two[0] = argrid_size[0] + 1;
        boats[i].arend_two[1] = argrid_size[1] + 1; //initialises while loops

        boats[i].ndir = rand_gen(1,4);


        switch(boats[i].ndir)
        {
        case 1:     //up
                while(boats[i].arend_two[1] > argrid_size[1] - 1)
                {
                    boats[i].nlength = rand_gen(2,4);   //arbitrarily set max/min boat sizes
                    boats[i].arend_two[0] = boats[i].arend_one[0];            //same x
                    boats[i].arend_two[1] = boats[i].arend_one[1] + boats[i].nlength - 1;    //bigger y
                }//loops until fits in grid
                for(short j = 0;j<boats[i].nlength;j++)
                    {
                        boats[i].boat_pos[1][j] = boats[i].arend_one[1] + j;
                        boats[i].boat_pos[0][j] = boats[i].arend_one[0];
                        //make array of boat positions
                    }
                break;
        case 2:     //right
                while(boats[i].arend_two[0] > argrid_size[0] - 1)
                {
                    boats[i].nlength = rand_gen(2,4);
                    boats[i].arend_two[1] = boats[i].arend_one[1];            //same y
                    boats[i].arend_two[0] = boats[i].arend_one[0] + boats[i].nlength - 1;    //bigger x

                }
                for(short j = 0;j<boats[i].nlength;j++)
                    {
                        boats[i].boat_pos[0][j] = boats[i].arend_one[0] + j;
                        boats[i].boat_pos[1][j] = boats[i].arend_one[1];
                        //make array of boat positions
                    }
                break;
        case 3:     //down
                while(boats[i].arend_two[1] < 0||boats[i].arend_two[1] > argrid_size[1] - 1)
                {
                    boats[i].nlength = rand_gen(2,4);
                    boats[i].arend_two[0] = boats[i].arend_one[0];            //same x
                    boats[i].arend_two[1] = boats[i].arend_one[1] - boats[i].nlength + 1;    //smaller y
                }
                for(short j = 0;j<boats[i].nlength;j++)
                    {
                        boats[i].boat_pos[1][j] = boats[i].arend_one[1] - j;
                        boats[i].boat_pos[0][j] = boats[i].arend_one[0];
                        //make array of boat positions
                    }
                break;
        case 4:     //left
                while(boats[i].arend_two[0] < 0||boats[i].arend_two[0] > argrid_size[0] - 1)
                {
                    boats[i].nlength = rand_gen(2,4);
                    boats[i].arend_two[1] = boats[i].arend_one[1];            //same y
                    boats[i].arend_two[0] = boats[i].arend_one[0] - boats[i].nlength + 1;     //smaller x
                }
                for(short j = 0;j<boats[i].nlength;j++)
                    {
                        boats[i].boat_pos[0][j] = boats[i].arend_one[0] + j;
                        boats[i].boat_pos[1][j] = boats[i].arend_one[1];
                        //make array of boat positions
                    }
                break;
        }
        //Calculates beginning and end points of the boat




    }
    //position boats in 10x10 matrix

    for(int i = 0;i < boat_number;i++)
        boats[i].boat_condition_set();
    //configure condition array for boats according to their length

    /*
    for(short k = 0;k<3;k++)
    {
     cout << "Length:\t";
     cout << boats[k].nlength << endl;

      for(short i = 0;i<2;i++)
    {
        cout << endl;
        for(short j = 0;j<4;j++)
            cout << boats[k].boat_pos[i][j] << endl;
    }
    }
    //walk through to help develop
    */


    while(james.ammo > 0)
    {
        james.receive_input();//is this the best way of doing this?
        if(hit_result(james,boats))
            cout << "Hit!" << endl;
        else
            cout << "Miss!" << endl;

        for(int i = 0;i<boat_number;i++)//sum over boats
            boats[i].afloat_check();

        sink_message(boats);
        if(victory_flag == true)
            return 0;

        if(james.ammo==20)
            cout << "20 ammo left!" << endl;
        else if(james.ammo==10)
            cout << "10 ammo left!!" << endl;
        else if(james.ammo==5)
            cout << "Only 5 ammo left!!!" << endl;
    }



    return 0;

}

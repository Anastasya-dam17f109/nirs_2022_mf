#include "pch.h"
#include "helix_curve.h"

// генерация точек кривой-спирали

void Helix_curve::get_points_for_curve(){
    Point cur_point = {0,0};
    bool stop_fl = false;
    bool up_fl = true;
    bool right_fl = true;
    int step_length = 1;
    switch(curve_type){
        case 0:{
            cur_point.x = n_order/2-1;
            cur_point.y = n_order/2;
            up_fl = true;
            right_fl = true;
        }
        break;
        case 1:{
            cur_point.x = n_order/2 - 1;
            cur_point.y = n_order/2 - 1;
            up_fl = false;
            right_fl = true;
        }
        break;
        case 2:{
            cur_point.x = n_order/2;
            cur_point.y = n_order/2 - 1;
            up_fl = false;
            right_fl = false;
        }
        break;
        case 3:{
            cur_point.x = n_order/2;
            cur_point.y = n_order/2;
            up_fl = true;
            right_fl = false;
        }
        break;
    }


    points.push_back(cur_point);
    while(!stop_fl){
        for(int i = 0; i< step_length; ++i){
            if(right_fl)
                cur_point.x++;
            else
                cur_point.x--;
            if((cur_point.x < 0)||(cur_point.x > (n_order-1))){
                stop_fl = true;
                break;
            }
            else
                points.push_back(cur_point);
        }

        if(!stop_fl){
            right_fl = !right_fl;
            for(int i = 0; i< step_length; ++i){
                if(up_fl)
                    cur_point.y--;
                else
                    cur_point.y++;
                if((cur_point.y < 0)||(cur_point.y > (n_order-1))){
                    stop_fl = true;
                    break;
                }
                else
                    points.push_back(cur_point);
            }
            up_fl = !up_fl;
            step_length++;
        }

    }

    if(reverse_flag)
        reverse_curve();
    //for(auto&i:points)
       // cout << "(" << i.x << "," << i.y << ")" <<endl;
}




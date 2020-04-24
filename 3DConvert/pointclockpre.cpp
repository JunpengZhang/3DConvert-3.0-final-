
#include "stdafx.h"
#include "pointclockpre.h"

//输入dx,dy为多边形顶点，rdx,rdy是按照顺序输出的多边形顶点
//对多边形顶点顺序进行判断，如果是逆时针则返回1，如果是顺时针则返回-1,输入不满足判断要求则返回0.
int ClockWisePro(int num, const double_vector & dx, const double_vector & dy, double_vector &rdx, double_vector &rdy)
{
    if ((num != dx.size()) || (num != dy.size()) || (num < 3))
    {
        return 0;
    }

    //计算有向面积，计算公式S = sum(X[i]*(Y[i+1] - Y[i-1])/2
    double area(0);
    area = dx[0]*(dy[1]- dy[num - 1])/2.0;
    for (int i = 1; i < num - 1; i++) 
    {
        area += dx[i]*(dy[i + 1] - dy[i - 1])/2.0;
    }

    area += dx[num - 1]*(dy[0] - dy[num - 2])/2.0;

    //面积S为负表明为顺时针，S非负为逆时针
    //顶点顺序为顺时针，直接返回原始数据
    if (area < 0.0)
    {
        rdx = dx;
        rdy = dy;
        return -1;
    }
    //顶点顺序为逆时针，转换为顺时针输出
    else
    {		
        for (int i  = 0; i < num; i++)
        {
            rdx.push_back(dx[num - 1 - i]);
            rdy.push_back(dy[num - 1 - i]);
        }
        return 1;
    }
}




//输入dx,dy为多边形顶点，rdx,rdy是按照逆序输出的多边形顶点
//对多边形顶点顺序进行判断，如果是逆时针则返回1，如果是顺时针则返回-1,输入不满足判断要求则返回0.
int antiClockWisePro(int num, const double_vector & dx, const double_vector & dy, double_vector &rdx, double_vector &rdy)
{
    if ((num != dx.size()) || (num != dy.size()) || (num < 3))
    {
        return 0;
    }

    //计算有向面积，计算公式S = sum(X[i]*(Y[i+1] - Y[i-1])/2
    double area(0);
    area = dx[0]*(dy[1]- dy[num - 1])/2.0;
    for (int i = 1; i < num - 1; i++) 
    {
        area += dx[i]*(dy[i + 1] - dy[i - 1])/2.0;
    }

    area += dx[num - 1]*(dy[0] - dy[num - 2])/2.0;

    //面积S为负表明为顺时针，S非负为逆时针
    //顶点顺序为逆时针，直接返回原始数据
    if (area > 0.0 || area == 0.0)
    {
        rdx = dx;
        rdy = dy;
        return 1;
    } 
    //顶点顺序为顺时针，转换为逆时针输出
    else
    {
        for (int i  = 0; i < num; i++)
        {
            rdx.push_back(dx[num - 1 - i]);
            rdy.push_back(dy[num - 1 - i]);
        }
        return -1;
    }	
}

#include "stdafx.h"
#include "pointclockpre.h"

//����dx,dyΪ����ζ��㣬rdx,rdy�ǰ���˳������Ķ���ζ���
//�Զ���ζ���˳������жϣ��������ʱ���򷵻�1�������˳ʱ���򷵻�-1,���벻�����ж�Ҫ���򷵻�0.
int ClockWisePro(int num, const double_vector & dx, const double_vector & dy, double_vector &rdx, double_vector &rdy)
{
    if ((num != dx.size()) || (num != dy.size()) || (num < 3))
    {
        return 0;
    }

    //����������������㹫ʽS = sum(X[i]*(Y[i+1] - Y[i-1])/2
    double area(0);
    area = dx[0]*(dy[1]- dy[num - 1])/2.0;
    for (int i = 1; i < num - 1; i++) 
    {
        area += dx[i]*(dy[i + 1] - dy[i - 1])/2.0;
    }

    area += dx[num - 1]*(dy[0] - dy[num - 2])/2.0;

    //���SΪ������Ϊ˳ʱ�룬S�Ǹ�Ϊ��ʱ��
    //����˳��Ϊ˳ʱ�룬ֱ�ӷ���ԭʼ����
    if (area < 0.0)
    {
        rdx = dx;
        rdy = dy;
        return -1;
    }
    //����˳��Ϊ��ʱ�룬ת��Ϊ˳ʱ�����
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




//����dx,dyΪ����ζ��㣬rdx,rdy�ǰ�����������Ķ���ζ���
//�Զ���ζ���˳������жϣ��������ʱ���򷵻�1�������˳ʱ���򷵻�-1,���벻�����ж�Ҫ���򷵻�0.
int antiClockWisePro(int num, const double_vector & dx, const double_vector & dy, double_vector &rdx, double_vector &rdy)
{
    if ((num != dx.size()) || (num != dy.size()) || (num < 3))
    {
        return 0;
    }

    //����������������㹫ʽS = sum(X[i]*(Y[i+1] - Y[i-1])/2
    double area(0);
    area = dx[0]*(dy[1]- dy[num - 1])/2.0;
    for (int i = 1; i < num - 1; i++) 
    {
        area += dx[i]*(dy[i + 1] - dy[i - 1])/2.0;
    }

    area += dx[num - 1]*(dy[0] - dy[num - 2])/2.0;

    //���SΪ������Ϊ˳ʱ�룬S�Ǹ�Ϊ��ʱ��
    //����˳��Ϊ��ʱ�룬ֱ�ӷ���ԭʼ����
    if (area > 0.0 || area == 0.0)
    {
        rdx = dx;
        rdy = dy;
        return 1;
    } 
    //����˳��Ϊ˳ʱ�룬ת��Ϊ��ʱ�����
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
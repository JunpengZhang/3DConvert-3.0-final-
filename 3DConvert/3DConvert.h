
// 3DConvert.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CMy3DConvertApp:
// �йش����ʵ�֣������ 3DConvert.cpp
//

class CMy3DConvertApp : public CWinApp
{
public:
	CMy3DConvertApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CMy3DConvertApp theApp;
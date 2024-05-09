#pragma once
#include "Scene.h"

namespace Tmpl8 {

class Surface;
class Game
{
public:
	void SetTarget( Surface* _Surface ) { screen = _Surface; }
	void Init();
	void Shutdown();
	void Tick( float _DT );
	void MouseUp( int _Button ) { /* implement if you want to detect mouse button presses */ }
	void MouseDown(int _Button, int _X, int _Y) {}
	void MouseMove(int _X, int _Y) { /* implement if you want to handle keys */ }
	void KeyUp( int _Key ) { /* implement if you want to handle keys */ }
	void KeyDown(int _Key);
	
	void loadObj(string obj);
	
	Surface* screen;
	Scene scene;
};

}; // namespace Tmpl8
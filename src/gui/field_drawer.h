#ifndef DUEL_FIELD_DRAWER_H_
#define DUEL_FIELD_DRAWER_H_

#include <memory>
#include <mutex>

#include "base/base.h"
#include "duel/game_state_observer.h"
#include "gui/box.h"
#include "gui/drawer.h"
#include "gui/screen.h"
#include "gui/unique_sdl_surface.h"
#include "gui/SDL_kanji.h"

class FieldRealtime;
class GameState;
class MainWindow;
class SDLCommentator;

// FieldDrawer draws the current puyo field etc.
class FieldDrawer : public GameStateObserver, public Drawer {
public:
    // Don't take ownership
    FieldDrawer();
    virtual ~FieldDrawer();

    virtual void onInit(Screen*) OVERRIDE;
    virtual void onUpdate(const GameState&) OVERRIDE;
    virtual void draw(Screen*) OVERRIDE;

private:
    void drawField(Screen*, int playerId, const FieldRealtime&);
    SDL_Rect toRect(PuyoColor);

    mutable std::mutex mu_;
    std::unique_ptr<GameState> gameState_;

    UniqueSDLSurface backgroundSurface_;
    UniqueSDLSurface puyoSurface_;
    UniqueSDLSurface ojamaSurface_;

    Kanji_Font* font_;
};

#endif

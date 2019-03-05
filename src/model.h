#ifndef MODEL_H
#define MODEL_H

extern const size_t K;
extern const size_t M;


double fullcond(const size_t comp, const double theta[K], size_t ngames,
                const size_t games[ngames][M], const size_t game_counts[ngames],
                const size_t win_counts[K]);

#endif

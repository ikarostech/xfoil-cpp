class BoundaryLayerUtil {
public:
    class EnvEnResult {
        public:
        double ax;
        double ax_hk;
        double ax_th;
        double ax_rt;
    };
    class AxResult {
        public:
        double ax;
        double ax_hk1;
        double ax_hk2;
        double ax_t1;
        double ax_t2;
        double ax_rt1;
        double ax_rt2;
        double ax_a1;
        double ax_a2;
    };
    static AxResult axset(double hk1, double thet1, double rt1, double a1, double hk2,
        double thet2, double rt2, double a2, double acrit);
    private:
    static EnvEnResult dampl(double hk, double thet, double rt);
};
class SlaterDeterminant{
  public:
    SlaterDeterminant();
    ~SlaterDeterminant();

    //SPS = Single Particle State
    SlaterDeterminant(int sps1, int sps2, int sps3, int sps4);
    std::vector<int> state_indices;
    int makeNormalOrder();
};

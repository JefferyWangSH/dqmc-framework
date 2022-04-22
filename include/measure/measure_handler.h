#ifndef MEASURE_HANDLER_H
#define MEASURE_HANDLER_H
#pragma once


/**
  *  
  *   
  */


#include "measure/observable_handler.h"


namespace Measure {


    // ----------------------------- Handler class Measure::MeasureHandler -----------------------------
    class MeasureHandler : public Observable::ObservableHandler {
        private:
            
            bool m_is_warmup{};
            bool m_is_equaltime{};
            bool m_is_dynamic{};

            int m_sweeps_warmup{};
            int m_bin_num{};
            int m_bin_size{};
            int m_sweeps_between_bins{};
        
        
        public:

            MeasureHandler() = default;

            const bool isEqualTime() const { return true; };
            const bool isDynamic() const { return true; };
            

            void initial() {};
    };


} // namespace Measure


#endif // MEASURE_HANDLER_H
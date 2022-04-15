#ifndef MEASURE_HANDLER_H
#define MEASURE_HANDLER_H
#pragma once


namespace Measure {

    // -------------------------- Handler class Measure::MeasureHandler -----------------------------
    class MeasureHandler {
        public:
            MeasureHandler() = default;

            const bool isEqualTime() const { return true; };
            const bool isDynamic() const { return true; };
            
    };


} // namespace Measure


#endif // MEASURE_HANDLER_H
//
//  aa_pipeline.hpp
//  aapipeline
//
//  Created by Cees Carels on Tuesday 23/10/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#ifndef ASTRO_ACCELERATE_PIPELINE_HPP
#define ASTRO_ACCELERATE_PIPELINE_HPP

#include <iostream>
#include <stdio.h>
#include <map>

#include "aa_compute.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_permitted_pipelines_2.hpp"

namespace astroaccelerate {

  /**
   * This class manages pipelines and their constituent modules, plans, and strategies,
   * as well as delegating the movement of host memory into and out of the pipeline,
   * and the movement of device memory into and out of the device.
   *
   * Still to do: Add a way to transfer ownership of the data between aa_pipeline objects.
   */

  template<typename T, typename U>
  class aa_pipeline {
  public:
    aa_pipeline(const aa_compute::pipeline &requested_pipeline, const aa_filterbank_metadata &filterbank_metadata, T const*const input_data, const aa_device_info::aa_card_info &card_info) : m_card_info(card_info), m_filterbank_metadata(filterbank_metadata), bound_with_raw_ptr(true), input_data_bound(true), data_on_device(false), pipeline_ready(false), ptr_data_in(input_data) {
        
      //Add requested pipeline modules
      for(auto i : requested_pipeline) {
	required_plans.insert(std::pair<aa_compute::modules, bool>(i, true));
	supplied_plans.insert(std::pair<aa_compute::modules, bool>(i, false));
      }
      m_all_strategy.reserve(requested_pipeline.size());
      m_requested_pipeline = requested_pipeline;
    }
    
    ~aa_pipeline() {
      pipeline_ready = false;
      if(data_on_device) {
	//There is still data on the device and the user has forgotten about it.
	//The device memory should be freed.
	std::cout << "ERROR: Data still on device." << std::endl;
      }
        
      if(!bound_with_raw_ptr && input_data_bound) {
	//There is still (host) input data bound to the pipeline and the user has forgotten about it.
	//This (host) memory should be freed.
	std::cout << "ERROR: Host data still bound to pipeline." << std::endl;
      }
    }
    
    //TODO: If any plan is offered to a bind method after a bind has already happened,
    //then all flags indicating readiness must be invalidated.
    bool bind(aa_ddtr_plan plan) {
      pipeline_ready = false;
        
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::dedispersion) != required_plans.end()) {
	m_ddtr_plan = plan;
            
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::dedispersion) = true;
            
	//ddtr_strategy needs to know if analysis will be required
	if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
	  aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, true);
	  if(ddtr_strategy.ready()) {
	    m_ddtr_strategy = std::move(ddtr_strategy);
	    m_all_strategy.push_back(&m_ddtr_strategy);
	  }
	  else {
	    return false;
	  }
	}
	else {
	  aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, false);
	  if(ddtr_strategy.ready()) {
	    m_ddtr_strategy = std::move(ddtr_strategy);
	    m_all_strategy.push_back(&m_ddtr_strategy);
	  }
	  else {
	    std::cout << "ddtr_strategy not ready" << std::endl;
	    return false;
	  }
	}
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }
    
    bool bind(aa_analysis_plan plan) {
      pipeline_ready = false;
        
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
	m_analysis_plan = plan;

	aa_analysis_strategy analysis_strategy(m_analysis_plan);
	if(analysis_strategy.ready()) {
	  m_analysis_strategy = std::move(analysis_strategy);
	  m_all_strategy.push_back(&m_analysis_strategy);
	}
	else {
	  return false;
	}
	
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::analysis) = true;
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }
    
    bool bind(aa_periodicity_plan plan) {
      pipeline_ready = false;
        
      //Does the pipeline actually need this plan?
      if(required_plans.find(aa_compute::modules::periodicity) != required_plans.end()) {
	m_periodicity_plan = plan;
	aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
	if(periodicity_strategy.ready()) {
	  m_periodicity_strategy = std::move(periodicity_strategy);
	  m_all_strategy.push_back(&m_periodicity_strategy);
	}
	else {
	  return false;
	}
	  
	//If the plan is valid then the supplied_plan becomes true
	supplied_plans.at(aa_compute::modules::periodicity) = true;
      }
      else {
	//The plan is not required, ignore.
	return false;
      }
        
      return true;
    }
    
    aa_ddtr_strategy ddtr_strategy() {
      //Does the pipeline actually need this strategy?
      if(required_plans.find(aa_compute::modules::dedispersion) != required_plans.end()) {
	//It does need this strategy.
	//Is it already computed?
	if(m_ddtr_strategy.ready()) {
	  //Return since it was already computed.
	  return m_ddtr_strategy;
	}
	else {
	  //ddtr_strategy was not yet computed, do it now.
	  //ddtr_strategy needs to know if analysis will be required
	  if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) { //analysis will be required
	    aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, true);
	    if(ddtr_strategy.ready()) {
	      m_ddtr_strategy = std::move(ddtr_strategy);
	      m_all_strategy.push_back(&m_ddtr_strategy);
	    }
	    else { //Tried to calculate ddtr_strategy with analysis enabled, but failed.
	      aa_ddtr_strategy empty_strategy;
	      return empty_strategy;
	    }
	  }
	  else { //analysis will not be required
	    aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_card_info.free_memory, false);
	    if(ddtr_strategy.ready()) {
	      m_ddtr_strategy = std::move(ddtr_strategy);
	      m_all_strategy.push_back(&m_ddtr_strategy);
	    }
	    else { //Tried to calculate ddtr_strategy with analysis disabled, but failed.
	      std::cout << "ddtr_strategy not ready" << std::endl;
	      aa_ddtr_strategy empty_strategy;
	      return empty_strategy;
	    }
	  }
	}
	
	return m_ddtr_strategy;
      }
      else {
	//The pipeline does not need this strategy
	aa_ddtr_strategy empty_strategy;
	return empty_strategy;
      }
    }

    aa_analysis_strategy analysis_strategy() {
      //Does the pipeline actually need this strategy? 
      if(required_plans.find(aa_compute::modules::analysis) != required_plans.end()) {
        //It does need this strategy.                                                                                                                                                                                        
        //Is it already computed?
        if(m_analysis_strategy.ready()) { //Return since it was already computed.
          return m_analysis_strategy;
        }
	else {
	  //analysis_strategy was not yet computed, do it now.
	  aa_analysis_strategy analysis_strategy(m_analysis_plan);
	  if(analysis_strategy.ready()) {
	    m_analysis_strategy = std::move(analysis_strategy);
	    m_all_strategy.push_back(&m_analysis_strategy);
	  }
	  else { //Tried to calculate analysis_strategy, but failed.
	    aa_analysis_strategy empty_strategy;
	    return empty_strategy;
	  }
	}
      }
      else {
	//The pipeline does not need this strategy
	aa_analysis_strategy empty_strategy;
	return empty_strategy;
      }
    }

    aa_periodicity_strategy periodicity_strategy() {
      //Does the pipeline actually need this strategy?
      if(required_plans.find(aa_compute::modules::periodicity) != required_plans.end()) {
	//It does need this strategy.
	//Is it already computed?
	if(m_periodicity_strategy.ready()) { //Return since it was already computed.
	  return m_periodicity_strategy;
	}
	else {
	  //periodicity_strategy was not yet computed, do it now.
	  aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
	  if(periodicity_strategy.ready()) {
	    m_periodicity_strategy = std::move(periodicity_strategy);
	    m_all_strategy.push_back(&m_periodicity_strategy);
	  }
	  else { //Tried to calculate periodicity strategy, but failed.
	    aa_periodicity_strategy empty_strategy;
	    return empty_strategy;
	  }
	}
      }
      else {
	//The pipeline does not need this strategy
	aa_periodicity_strategy empty_strategy;
	return empty_strategy;
	
      }
    }
    
    bool transfer_data_to_device() {
      pipeline_ready = false;
      if(input_data_bound) {
	data_on_device = true;
      }
        
      return true;
    }
    
    bool transfer_data_to_host(std::vector<U> &data) {
      if(data_on_device) {
	data = std::move(data_out);
	data_on_device = false;
	pipeline_ready = false;
	return true;
      }
        
      return false;
    }
    
    bool transfer_data_to_host(U *&data) {
      if(data_on_device) {
	data = ptr_data_out;
	data_on_device = false;
	pipeline_ready = false;
	return true;
      }
        
      return false;
    }
    
    bool unbind_data() {
      /**
       * If the data is managed, then either it was either moved out via the transfer,
       * or it will be freed at de-allocation if the user forgot to transfer.
       *
       * If the data is unmanaged, then the data unbind call does not apply.
       * If the data came via a raw pointer, then the unbind call does not apply.
       *
       */
      pipeline_ready = false;
      return true;
    }
    
    bool ready() {
      pipeline_ready = false;
        
      if(!data_on_device) {
	return false;
      }
        
      for(auto const& i : supplied_plans) {
	if(i.second == false) {
	  return false;
	}
      }
        
      //Do any last checks on the plans as a whole
      pipeline_ready = true;
        
      return true;
    }
    
    bool run() {
      if(pipeline_ready && data_on_device) {
	//A function callback can run the pipeline from elsewhere
	for(auto strategy : m_all_strategy) {
	  if(strategy->setup()) {
	    //Memory allocations and setup happened successfully
	  }
	  else {
	    //Setup failed.
	  }
	}

	//Choose a specialisation of zero_dm for now
	aa_permitted_pipelines_1<aa_compute::modules::zero_dm, false> runner(m_ddtr_strategy, ptr_data_in);
	if(runner.setup()) {
	  std::vector<float> out;
	  int chunk_idx = 0;
	  std::vector<int> range_samples;
	  while(runner.next(out, chunk_idx, range_samples)) {
	    std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
	  }
	}
            
	return true;
      }
      else {
	return false;
      }
    }
    
    bool handoff(aa_pipeline &next_pipeline) {
      /**
       * Handoff control over the data to the next pipeline.
       *
       * The user must always retrieve the data from the last pipeline in the chain.
       */
      return true;
    }
    
  private:
    std::map<aa_compute::modules, bool> required_plans; //Plans required to configure the pipeline
    std::map<aa_compute::modules, bool> supplied_plans; //Plans supplied by the user
    std::vector<aa_strategy*>   m_all_strategy;
    aa_compute::pipeline        m_requested_pipeline;
    aa_device_info::aa_card_info m_card_info;
    
    aa_filterbank_metadata      m_filterbank_metadata;
    
    aa_ddtr_plan                m_ddtr_plan;
    aa_ddtr_strategy            m_ddtr_strategy;
    
    aa_analysis_plan            m_analysis_plan;
    aa_analysis_strategy        m_analysis_strategy;
    
    aa_periodicity_plan         m_periodicity_plan;
    aa_periodicity_strategy     m_periodicity_strategy;
    
    bool bound_with_raw_ptr;
    bool input_data_bound;
    bool data_on_device;
    bool pipeline_ready;
    
    std::vector<T>              data_in;
    std::vector<U>              data_out;
    T const*const               ptr_data_in;
    U*                          ptr_data_out;
  };

} //namespace astroaccelerate
  
#endif /* ASTRO_ACCELERATE_PIPELINE_HPP */

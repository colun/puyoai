#include "core/server/connector/human_connector.h"

#include <glog/logging.h>

#include "core/frame_request.h"
#include "core/frame_response.h"

using namespace std;

void HumanConnector::send(const FrameRequest& req)
{
    writeString(req.toString());
}

void HumanConnector::writeString(const string& message)
{
    LOG(INFO) << message;
}

bool HumanConnector::receive(FrameResponse* response)
{
    *response = FrameResponse();

    lock_guard<mutex> lock(mu_);
    response->keySet = currentKeySet_;
    return true;
}

void HumanConnector::setClosed(bool)
{
    CHECK(false) << "HumanConnector does not have closed flag.";
}

int HumanConnector::readerFd() const
{
    CHECK(false) << "HumanConnector does not have reader file descriptor.";
    return -1;
}

void HumanConnector::setKeySet(const KeySet& keySet)
{
    lock_guard<mutex> lock(mu_);
    currentKeySet_ = keySet;
}
